#' WriteH5AD_Spatial
#' @description
#' Write a spatial Seurat object into h5ad AnnData format, preserving spatial
#' coordinates, tissue images, and scale factors.
#'
#' @details
#' Supports Seurat objects with VisiumV1 or VisiumV2 image objects.
#' Spatial data is written to standard AnnData locations:
#' \itemize{
#'   \item \code{obsm["spatial"]}: Tissue coordinates (imagecol, imagerow),
#'     following the Scanpy convention of (x, y)
#'   \item \code{uns["spatial"][library_id]["images"]}: H&E image as numpy array
#'   \item \code{uns["spatial"][library_id]["scalefactors"]}: Scale factor dictionary
#'   \item \code{obs}: Includes in_tissue, array_row, array_col columns
#' }
#'
#' @param seurat_object The Seurat object with spatial data to be converted
#' @param env_path A path used to execute python environment.
#' @param output_path The path to store the converted AnnData object
#' @param assay The name of the assay to use for expression data. Default "Spatial".
#' @param image_name The name of the image/slice in the Seurat object. If NULL,
#'   the first available image is used.
#' @param library_id The library ID to use in \code{uns["spatial"]}. Default "library1".
#' @param image_resolution Which image resolution to write: "hires", "lowres", or "both". Default "lowres",
#'   matching the image normally stored in a Seurat Visium object.
#' @param hires_image_path Optional path (character) to an external hires tissue
#'   image PNG file (e.g. Space Ranger's \code{tissue_hires_image.png}).
#'   Default \code{NULL}. Behaviour:
#'   \itemize{
#'     \item If supplied while \code{image_resolution} is omitted, both hires
#'       and the in-object lowres image are written when available.
#'     \item If \code{NULL}: the hires image is taken from the Seurat object's
#'       \code{@image} slot when that slot contains a hires-sized image
#'       (detected by dimension; typically >=1500 px on the long edge).
#'       Otherwise the slot is treated as lowres.
#'     \item If a character path: the file at that path is read with
#'       \code{png::readPNG} and used as the hires image, overriding whatever
#'       is in the Seurat object.
#'   }
#'
#' @examples
#' \dontrun{
#'   WriteH5AD_Spatial(spatial_obj, "path/to/python/env", "spatial.h5ad")
#'   WriteH5AD_Spatial(spatial_obj, "path/to/python/env", "spatial.h5ad",
#'                     image_resolution = "both",
#'                     hires_image_path = "spatial/tissue_hires_image.png")
#' }
#'
#' @importFrom reticulate use_python import py_to_r py_call
#' @importFrom Seurat Embeddings Images GetAssayData
#' @importFrom methods as is slot
#' @importFrom Matrix t
#' @return No return value. The AnnData object is written to \code{output_path}.
#' @export

WriteH5AD_Spatial = function(seurat_object, env_path, output_path,
                              assay = "Spatial", image_name = NULL,
                              library_id = "library1",
                              image_resolution = "lowres",
                              hires_image_path = NULL) {

  if (!missing(env_path)) {
    use_python(env_path, required = TRUE)
  }

  if (!is.null(hires_image_path) && missing(image_resolution)) {
    image_resolution = "both"
  }

  if (!(image_resolution %in% c("hires", "lowres", "both"))) {
    stop("image_resolution must be one of 'hires', 'lowres', 'both'.")
  }
  if (!is.character(library_id) || length(library_id) != 1 ||
      is.na(library_id) || !nzchar(library_id)) {
    stop("library_id must be a single non-empty string.")
  }

  load_image_file = function(path) {
    if (!grepl("\\.png$", path, ignore.case = TRUE)) {
      stop("hires_image_path must be a .png file (got '", path,
           "'). Space Ranger writes tissue_hires_image.png by default.")
    }
    if (!requireNamespace("png", quietly = TRUE)) {
      stop("Package 'png' is required to read '", path,
           "'. Install it with install.packages('png').")
    }
    png::readPNG(path)
  }

  # External hires image (if user passed a file path) — wins over whatever is
  # stored in the Seurat object. NULL means: try to recover hires from the
  # @image slot via dimension-based detection further below.
  external_hires_img = NULL
  if (!is.null(hires_image_path)) {
    if (!is.character(hires_image_path) || length(hires_image_path) != 1) {
      stop("hires_image_path must be NULL or a single file path string.")
    }
    if (!file.exists(hires_image_path)) {
      stop("hires_image_path does not exist: ", hires_image_path)
    }
    external_hires_img = load_image_file(hires_image_path)
  }

  # Heuristic for classifying the in-object image. Visium hires images are
  # typically ~2000 px on the long edge; lowres are ~600 px. Anything >=1500
  # px is treated as hires.
  classify_slot_image = function(img) {
    if (is.null(img) || !is.array(img) || length(dim(img)) != 3) return(NA_character_)
    if (max(dim(img)[1:2]) >= 1500) "hires" else "lowres"
  }

  tryCatch({
    anndata = reticulate::import("anndata", delay_load = FALSE)
    np = reticulate::import("numpy", delay_load = FALSE)
    pd = reticulate::import("pandas", delay_load = FALSE)
  }, error = function(e) {
    stop("Required Python packages (anndata, numpy, pandas) not found: ", e$message)
  })

  # Define a Python helper that performs all spatial assignments directly on
  # the AnnData object. Doing this in Python avoids reticulate's auto-conversion
  # of `adata.uns` / `adata.obsm` to detached R lists, which is why R-side
  # `[<-` and even `operator.setitem` calls failed to persist.
  reticulate::py_run_string("
def _gauguin_set_spatial(adata, coords, library_id,
                         hires_img=None, lowres_img=None,
                         scalefactors=None, source_image_path=''):
    import numpy as _np
    adata.obsm['spatial'] = _np.asarray(coords)
    lib = {'metadata': {'source_image_path': str(source_image_path)}}
    images = {}
    if hires_img is not None:
        images['hires'] = _np.asarray(hires_img)
    if lowres_img is not None:
        images['lowres'] = _np.asarray(lowres_img)
    if images:
        lib['images'] = images
    if scalefactors is not None:
        lib['scalefactors'] = {str(k): v for k, v in dict(scalefactors).items()}
    adata.uns['spatial'] = {str(library_id): lib}

def _gauguin_set_obsm(adata, key, value):
    import numpy as _np
    adata.obsm[str(key)] = _np.asarray(value)
")
  py_set_spatial = reticulate::py$`_gauguin_set_spatial`
  py_set_obsm    = reticulate::py$`_gauguin_set_obsm`

  # Detect available assay: try user-specified, fallback to common spatial assay names
  available_assays = names(seurat_object@assays)
  if (!(assay %in% available_assays)) {
    requested_assay = assay
    spatial_candidates = c("Spatial", "SCT", "RNA")
    matched = spatial_candidates[spatial_candidates %in% available_assays]
    if (length(matched) > 0) {
      assay = matched[1]
      message("Assay '", requested_assay, "' not found, using '", assay, "' instead.")
    } else {
      assay = available_assays[1]
      message("Using first available assay: '", assay, "'")
    }
  }

  # Extract counts matrix
  counts_matrix = Seurat::GetAssayData(seurat_object, assay = assay, layer = "counts")
  if (!inherits(counts_matrix, "matrix") && !inherits(counts_matrix, "dgCMatrix")) {
    counts_matrix = as.matrix(counts_matrix)
  }

  gene_names = rownames(counts_matrix)
  cell_names = colnames(counts_matrix)
  counts_matrix = Matrix::t(counts_matrix)
  counts_matrix = as(counts_matrix, "dgCMatrix")

  # Metadata
  meta_data = seurat_object@meta.data
  if (is.null(rownames(meta_data)) || !all(cell_names %in% rownames(meta_data))) {
    stop("Seurat metadata row names do not contain all assay cell barcodes.")
  }
  meta_data = meta_data[cell_names, , drop = FALSE]
  meta_data$barcode = cell_names

  create_anndata = function(counts, metadata, cells) {
    metadata = metadata[cells, , drop = FALSE]
    adata = anndata$AnnData(X = counts, obs = pd$DataFrame(metadata))
    adata$obs_names = np$array(cells)
    adata$var_names = np$array(gene_names)
    adata
  }

  write_reductions = function(adata, cells) {
    if (length(seurat_object@reductions) == 0) return(invisible(NULL))

    for (reduction_name in names(seurat_object@reductions)) {
      embeddings = Embeddings(seurat_object, reduction = reduction_name)
      if (is.null(rownames(embeddings)) || !all(cells %in% rownames(embeddings))) {
        warning("Skipping reduction '", reduction_name,
                "' because its cell names do not match the exported cells.")
        next
      }
      embeddings = embeddings[cells, , drop = FALSE]
      obsm_key = if (startsWith(reduction_name, "X_")) {
        reduction_name
      } else {
        paste0("X_", reduction_name)
      }
      py_set_obsm(adata, obsm_key, unname(embeddings))
    }
    invisible(NULL)
  }

  # --- Spatial data migration ---
  image_names = Seurat::Images(seurat_object)

  if (length(image_names) == 0) {
    warning("No spatial images found in Seurat object. Writing non-spatial H5AD.")
    adata = create_anndata(counts_matrix, meta_data, cell_names)
    write_reductions(adata, cell_names)
    adata$write_h5ad(output_path)
    message("Saved (non-spatial) at ", output_path)
    return(invisible(NULL))
  }

  # Select image
  if (is.null(image_name)) {
    image_name = image_names[1]
  }
  if (!(image_name %in% image_names)) {
    stop("Image '", image_name, "' not found. Available: ", paste(image_names, collapse = ", "))
  }

  image_obj = seurat_object[[image_name]]
  message("Migrating spatial data from image '", image_name, "'...")

  # --- Extract coordinates based on image class ---
  # Note: obsm["spatial"] in scanpy/squidpy is (x, y) = (imagecol, imagerow).
  # Seurat's `@image` slot stores the *lowres* tissue image, not the hires one.
  spatial_coords = NULL
  tissue_vec = NULL
  array_row = NULL
  array_col = NULL
  lowres_img_array = NULL
  sf_list = NULL

  if (is(image_obj, "VisiumV1")) {
    # VisiumV1: coordinates in @coordinates slot
    coords_df = slot(image_obj, "coordinates")

    spatial_coords = as.matrix(coords_df[, c("imagecol", "imagerow")])

    if ("tissue" %in% colnames(coords_df)) {
      tissue_vec = as.integer(coords_df[["tissue"]])
    }
    if ("row" %in% colnames(coords_df)) {
      array_row = as.integer(coords_df[["row"]])
    }
    if ("col" %in% colnames(coords_df)) {
      array_col = as.integer(coords_df[["col"]])
    }

    # Image (lowres — that's what Seurat keeps)
    lowres_img_array = slot(image_obj, "image")

    # Scale factors
    sf_obj = slot(image_obj, "scale.factors")
    sf_list = list(
      spot_diameter_fullres = sf_obj$spot,
      fiducial_diameter_fullres = sf_obj$fiducial,
      tissue_hires_scalef = sf_obj$hires,
      tissue_lowres_scalef = sf_obj$lowres
    )

  } else if (is(image_obj, "VisiumV2")) {
    # VisiumV2: inherits FOV; GetTissueCoordinates returns columns (x, y, cell)
    centroids = SeuratObject::GetTissueCoordinates(image_obj)
    if (all(c("x", "y") %in% colnames(centroids))) {
      spatial_coords = as.matrix(centroids[, c("x", "y")])
    } else {
      spatial_coords = as.matrix(centroids[, c(1, 2)])
    }
    cell_column = intersect(c("cell", "barcode"), colnames(centroids))
    if (length(cell_column) > 0) {
      rownames(spatial_coords) = as.character(centroids[[cell_column[1]]])
    }

    # Image (lowres) and scale factors from VisiumV2 slots
    lowres_img_array = tryCatch(slot(image_obj, "image"), error = function(e) NULL)
    sf_obj = tryCatch(slot(image_obj, "scale.factors"), error = function(e) NULL)
    if (!is.null(sf_obj)) {
      sf_list = list(
        spot_diameter_fullres = sf_obj$spot,
        fiducial_diameter_fullres = sf_obj$fiducial,
        tissue_hires_scalef = sf_obj$hires,
        tissue_lowres_scalef = sf_obj$lowres
      )
    }

  } else if (is(image_obj, "FOV")) {
    # Generic FOV (Xenium, Vizgen, etc.): centroids from GetTissueCoordinates
    centroids = SeuratObject::GetTissueCoordinates(image_obj)
    if (all(c("x", "y") %in% colnames(centroids))) {
      spatial_coords = as.matrix(centroids[, c("x", "y")])
    } else {
      spatial_coords = as.matrix(centroids[, c(1, 2)])
    }
    cell_column = intersect(c("cell", "barcode"), colnames(centroids))
    if (length(cell_column) > 0) {
      rownames(spatial_coords) = as.character(centroids[[cell_column[1]]])
    }

  } else if (is(image_obj, "SlideSeq")) {
    coords_df = slot(image_obj, "coordinates")
    spatial_coords = as.matrix(coords_df[, 1:2])
  }

  if (is.null(spatial_coords)) {
    warning("Could not extract spatial coordinates from image object of class '",
            class(image_obj)[1], "'. Writing non-spatial H5AD.")
    adata = create_anndata(counts_matrix, meta_data, cell_names)
    write_reductions(adata, cell_names)
    adata$write_h5ad(output_path)
    message("Saved (non-spatial) at ", output_path)
    return(invisible(NULL))
  }

  if (!is.null(sf_list)) {
    valid_scale = vapply(
      sf_list,
      function(value) {
        is.numeric(value) && length(value) == 1 &&
          is.finite(value) && value > 0
      },
      logical(1)
    )
    if (any(!valid_scale)) {
      warning(
        "Omitting invalid spatial scale factors: ",
        paste(names(sf_list)[!valid_scale], collapse = ", ")
      )
      sf_list = sf_list[valid_scale]
    }
    if (length(sf_list) == 0) sf_list = NULL
  }

  # Name the index-aligned obs vectors with the same row order as spatial_coords
  # so we can subset everything together when cells don't fully overlap.
  coord_rownames = rownames(spatial_coords)
  if (!is.null(coord_rownames) && anyDuplicated(coord_rownames)) {
    stop("Spatial coordinate cell barcodes must be unique.")
  }
  name_vec = function(v) {
    if (!is.null(v) && !is.null(coord_rownames) && length(v) == length(coord_rownames)) {
      names(v) = coord_rownames
    }
    v
  }
  tissue_vec = name_vec(tissue_vec)
  array_row  = name_vec(array_row)
  array_col  = name_vec(array_col)

  # Ensure coordinates align with cells in the counts matrix
  common_cells = cell_names[cell_names %in% coord_rownames]
  if (length(common_cells) == 0) {
    # Positional matching is only safe when coordinates have no meaningful
    # row names. Never relabel an equally sized but conflicting barcode set.
    default_coord_names = is.null(coord_rownames) || identical(
      coord_rownames,
      as.character(seq_len(nrow(spatial_coords)))
    )
    if (default_coord_names && nrow(spatial_coords) == length(cell_names)) {
      rownames(spatial_coords) = cell_names
      coord_rownames = cell_names
      if (!is.null(tissue_vec)) names(tissue_vec) = cell_names
      if (!is.null(array_row))  names(array_row)  = cell_names
      if (!is.null(array_col))  names(array_col)  = cell_names
      common_cells = cell_names
    } else {
      warning("Spatial coordinates do not match cell barcodes. Skipping spatial data.")
      adata = create_anndata(counts_matrix, meta_data, cell_names)
      write_reductions(adata, cell_names)
      adata$write_h5ad(output_path)
      message("Saved (non-spatial) at ", output_path)
      return(invisible(NULL))
    }
  }

  # If overlap is partial, drop non-overlapping cells from counts and metadata
  # so adata, meta_data and spatial_coords all line up on `common_cells`.
  if (length(common_cells) < length(cell_names)) {
    message("Spatial coords cover ", length(common_cells), "/", length(cell_names),
            " cells; dropping unmatched cells from H5AD.")
    counts_matrix = counts_matrix[common_cells, , drop = FALSE]
    meta_data     = meta_data[common_cells, , drop = FALSE]
    cell_names    = common_cells
  }

  # Reorder coordinates and obs vectors to match adata's cell order
  spatial_coords = spatial_coords[cell_names, , drop = FALSE]
  if (!is.null(tissue_vec)) tissue_vec = tissue_vec[cell_names]
  if (!is.null(array_row))  array_row  = array_row[cell_names]
  if (!is.null(array_col))  array_col  = array_col[cell_names]

  # Inject tissue-position obs columns BEFORE constructing AnnData so the
  # pandas DataFrame ends up with the columns (DataFrame proxy assignment
  # from R is unreliable).
  if (!is.null(tissue_vec) && !("in_tissue" %in% colnames(meta_data))) {
    meta_data[["in_tissue"]] = as.integer(tissue_vec)
  }
  if (!is.null(array_row) && !("array_row" %in% colnames(meta_data))) {
    meta_data[["array_row"]] = as.integer(array_row)
  }
  if (!is.null(array_col) && !("array_col" %in% colnames(meta_data))) {
    meta_data[["array_col"]] = as.integer(array_col)
  }
  if (!("library_id" %in% colnames(meta_data))) {
    meta_data[["library_id"]] = library_id
  }

  # Create AnnData now that meta_data is final
  adata = create_anndata(counts_matrix, meta_data, cell_names)

  # --- Dimensional reductions ---
  write_reductions(adata, cell_names)

  # Decide which image resolutions to include.
  # Priority for hires: external file (hires_image_path) > in-object @image
  # slot when its dimensions look like hires. lowres comes from the @image
  # slot when its dimensions look like lowres.
  is_valid_img = function(x) !is.null(x) && is.array(x) && length(dim(x)) == 3

  slot_role = classify_slot_image(lowres_img_array)  # "hires" | "lowres" | NA

  hires_from_slot  = if (identical(slot_role, "hires"))  lowres_img_array else NULL
  lowres_from_slot = if (identical(slot_role, "lowres")) lowres_img_array else NULL

  hires_arg = NULL
  if (image_resolution %in% c("hires", "both")) {
    if (is_valid_img(external_hires_img)) {
      hires_arg = external_hires_img
    } else if (is_valid_img(hires_from_slot)) {
      hires_arg = hires_from_slot
      message("  Using hires image found in Seurat @image slot (",
              paste(dim(hires_from_slot)[1:2], collapse = "x"), " px).")
    } else {
      warning("image_resolution requested '", image_resolution,
              "' but no hires image is available (hires_image_path is NULL ",
              "and the Seurat @image slot does not look like hires). ",
              "Hires will be omitted; spot overlay should still work via lowres.")
    }
  }

  lowres_arg = NULL
  if (image_resolution %in% c("lowres", "both") && is_valid_img(lowres_from_slot)) {
    lowres_arg = lowres_from_slot
  }

  # All spatial writes happen inside the Python helper to ensure they mutate
  # the live AnnData object (reticulate may auto-convert adata.uns / adata.obsm
  # to detached R lists when assigned from the R side).
  py_set_spatial(
    adata        = adata,
    coords       = unname(spatial_coords),
    library_id   = library_id,
    hires_img    = hires_arg,
    lowres_img   = lowres_arg,
    scalefactors = sf_list,
    source_image_path = if (is.null(hires_image_path)) "" else normalizePath(
      hires_image_path, winslash = "/", mustWork = TRUE
    )
  )
  message("  Written obsm['spatial'] and uns['spatial'][", library_id, "]")

  # Write output
  adata$write_h5ad(output_path)
  message("Spatial H5AD saved at ", output_path)
  return(invisible(NULL))
}
