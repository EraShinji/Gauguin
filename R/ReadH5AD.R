#' ReadH5AD
#' @description
#' Read H5AD AnnData Object, and convert to Seurat Object.
#' Supports spatial data (Visium) by automatically detecting \code{obsm["spatial"]}
#' and \code{uns["spatial"]} and constructing a VisiumV1 image object.
#'
#' @param h5ad_path A path which stores H5AD file.
#' @param env_path A path used to execute python environment.
#' @returns seurat_obj A Seurat Object
#' @examples
#' \dontrun{
#'
#' h5ad_file = system.file("extdata", "sce.h5ad", package = "Gauguin")
#' ReadH5AD(h5ad_file, env_path = "path/to/your/python/env")
#' }
#' @importFrom reticulate use_python import py_to_r py_is_null_xptr py_call import_builtins
#' @importFrom Seurat CreateSeuratObject SetAssayData CreateDimReducObject
#' @importFrom methods as new
#' @importFrom Matrix t
#' @export
ReadH5AD = function(h5ad_path, env_path) {

  if (!missing(env_path)) {
    use_python(env_path, required = TRUE)
  }

  anndata = import("anndata", convert = FALSE)
  adata = anndata$read_h5ad(h5ad_path)

  obs_names = as.character(py_to_r(adata$obs_names$to_list()))
  var_names = as.character(py_to_r(adata$var_names$to_list()))
  builtins <- import_builtins()
  layers_list <- as.character(py_to_r(builtins$list(adata$layers$keys())))

  if ("counts" %in% layers_list) {
    message("Using counts layer...")
    counts_mtx = adata$layers["counts"]
  } else if (!py_is_null_xptr(adata$raw)) {
    message("Using adata.raw.X...")
    counts_mtx = adata$raw$X
    var_names = as.character(py_to_r(adata$raw$var_names$to_list()))
  } else {
    message("Using adata.X as counts...")
    counts_mtx = adata$X
  }

  counts_mtx = py_to_r(counts_mtx)
  if (!inherits(counts_mtx, "CsparseMatrix")) {
    counts_mtx = as(counts_mtx, "CsparseMatrix")
  }
  counts_mtx = Matrix::t(counts_mtx)

  rownames(counts_mtx) = var_names
  colnames(counts_mtx) = obs_names

  meta_data = py_to_r(adata$obs)
  rownames(meta_data) = obs_names

  seurat_obj = CreateSeuratObject(counts = counts_mtx, meta.data = meta_data, assay = "RNA")


  if ("data" %in% layers_list) {
    message("Migrating 'data' layer...")
    data_mtx = py_to_r(adata$layers["data"])
    data_mtx = Matrix::t(data_mtx)

    rownames(data_mtx) = var_names
    colnames(data_mtx) = obs_names

    seurat_obj = SetAssayData(seurat_obj, layer = "data", new.data = as(data_mtx, "CsparseMatrix"))

  } else if (!py_is_null_xptr(adata$X)) {
    message(" cannot find `data` layer in adata.layers, try to use adata.X if normalized instead")

    n_obs_r = as.integer(py_to_r(adata$n_obs))
    n_vars_r = as.integer(py_to_r(adata$n_vars))

    test_vals = py_to_r(adata$X[0:min(5, n_obs_r), 0:min(100, n_vars_r)])
    is_float = any(as.numeric(test_vals) %% 1 != 0)

    if (is_float) {

      message("adata.X contains float values, migrating to 'data' layer...")
      data_mtx = Matrix::t(py_to_r(adata$X))
      rownames(data_mtx) = var_names
      colnames(data_mtx) = obs_names
      seurat_obj = SetAssayData(seurat_obj, layer = "data", new.data = as(data_mtx, "CsparseMatrix"))
    } else {
      message("adata.X appears to be integer counts, skipping 'data' layer migration.")
    }

  }

  n_obs_r = as.integer(py_to_r(adata$n_obs))
  n_vars_r = as.integer(py_to_r(adata$n_vars))
  sample_X = py_to_r(adata$X[0:min(100, n_obs_r-1), 0:min(100, n_vars_r-1)])


  is_scaled = any(as.numeric(sample_X) < -0.0001)
  if ("scaled" %in% layers_list) {
    scaled_mtx = Matrix::t(py_to_r(adata$layers["scaled"]))
    rownames(scaled_mtx) = as.character(py_to_r(adata$var_names$to_list()))
    colnames(scaled_mtx) = obs_names
    seurat_obj = SetAssayData(seurat_obj, slot = "scale.data", new.data = as.matrix(scaled_mtx))
  } else if(is_scaled){
    message("adata.X might be scaled, please consider use adata.X as layer `scale.data`")
    scaled_mtx = Matrix::t(py_to_r(adata$X))
  }else{
    message("No scaled.data be yieled in anndata object,nor any scaled matrix be generated")
  }

  obsm_dict = py_call(adata$obsm$as_dict)

  for (key in names(obsm_dict)) {
    if (key == "spatial") next
    embedding = obsm_dict[[key]]
    embedding = py_to_r(embedding)
    rownames(embedding) = colnames(seurat_obj)
    colnames(embedding) = paste0(key, "_", 1:ncol(embedding))
    fixed_key = gsub("^X_", "", key)
    fixed_key = paste0(fixed_key, "")
    seurat_obj[[fixed_key]] = CreateDimReducObject(embeddings = embedding,
                                                    key = fixed_key, assay = "RNA")
  }

  # --- Spatial data migration ---
  has_spatial_obsm = "spatial" %in% names(obsm_dict)
  uns_dict = tryCatch(py_to_r(py_call(adata$uns$as_dict)), error = function(e) list())
  has_spatial_uns = "spatial" %in% names(uns_dict)

  if (has_spatial_obsm) {
    message("Detected spatial data in obsm['spatial'], migrating...")

    spatial_coords = py_to_r(obsm_dict[["spatial"]])
    rownames(spatial_coords) = obs_names

    # Build coordinates data.frame
    # AnnData spatial coords: column 1 = imagerow (pxl_row), column 2 = imagecol (pxl_col)
    coordinates = data.frame(
      imagerow = spatial_coords[, 1],
      imagecol = spatial_coords[, 2],
      row.names = obs_names
    )

    # Add tissue, array_row, array_col from obs if available
    obs_df = py_to_r(adata$obs)
    if ("in_tissue" %in% colnames(obs_df)) {
      coordinates$tissue = as.integer(obs_df[["in_tissue"]])
    } else {
      coordinates$tissue = 1L
    }
    if ("array_row" %in% colnames(obs_df)) {
      coordinates$row = as.integer(obs_df[["array_row"]])
    } else {
      coordinates$row = seq_len(nrow(coordinates))
    }
    if ("array_col" %in% colnames(obs_df)) {
      coordinates$col = as.integer(obs_df[["array_col"]])
    } else {
      coordinates$col = seq_len(nrow(coordinates))
    }

    # Reorder columns to match Seurat convention
    coordinates = coordinates[, c("tissue", "row", "col", "imagerow", "imagecol")]

    # Extract image and scale factors from uns["spatial"] if available
    spatial_image = array(0, dim = c(1, 1, 3))
    sf_spot = 1
    sf_fiducial = 1
    sf_hires = 1
    sf_lowres = 1

    if (has_spatial_uns) {
      spatial_uns = uns_dict[["spatial"]]
      library_id = names(spatial_uns)[1]

      if (!is.null(library_id)) {
        lib_data = spatial_uns[[library_id]]

        # Extract images
        if ("images" %in% names(lib_data)) {
          images_data = lib_data[["images"]]
          if ("hires" %in% names(images_data)) {
            spatial_image = images_data[["hires"]]
            message("  Loaded hires image from uns['spatial']")
          } else if ("lowres" %in% names(images_data)) {
            spatial_image = images_data[["lowres"]]
            message("  Loaded lowres image from uns['spatial']")
          }
          # Ensure image is a 3D array
          if (!is.array(spatial_image)) {
            spatial_image = as.array(spatial_image)
          }
        }

        # Extract scale factors
        if ("scalefactors" %in% names(lib_data)) {
          sf_data = lib_data[["scalefactors"]]
          if ("spot_diameter_fullres" %in% names(sf_data))
            sf_spot = as.numeric(sf_data[["spot_diameter_fullres"]])
          if ("fiducial_diameter_fullres" %in% names(sf_data))
            sf_fiducial = as.numeric(sf_data[["fiducial_diameter_fullres"]])
          if ("tissue_hires_scalef" %in% names(sf_data))
            sf_hires = as.numeric(sf_data[["tissue_hires_scalef"]])
          if ("tissue_lowres_scalef" %in% names(sf_data))
            sf_lowres = as.numeric(sf_data[["tissue_lowres_scalef"]])
          message("  Loaded scale factors from uns['spatial']")
        }
      }
    }

    # Create scalefactors object
    scale_factors = structure(
      list(spot = sf_spot, fiducial = sf_fiducial, hires = sf_hires, lowres = sf_lowres),
      class = "scalefactors"
    )

    # Create VisiumV1 object
    visium_image = new(
      Class = "VisiumV1",
      assay = "RNA",
      key = "slice1_",
      image = spatial_image,
      scale.factors = scale_factors,
      coordinates = coordinates,
      spot.radius = sf_spot / 2
    )

    # Attach image to Seurat object
    seurat_obj[["slice1"]] = visium_image
    message("Spatial data successfully attached as VisiumV1 image 'slice1'")
  }

  return(seurat_obj)
}
