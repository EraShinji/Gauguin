test_that("WriteH5AD_Spatial preserves barcodes, coordinate order, images, scales, and reductions", {
  skip_if(!reticulate::py_module_available("anndata"), "anndata is not available")

  cells <- paste0("cell", 1:4)
  counts <- Matrix::Matrix(
    matrix(
      1:12,
      nrow = 3,
      dimnames = list(paste0("gene", 1:3), cells)
    ),
    sparse = TRUE
  )
  object <- Seurat::CreateSeuratObject(counts = counts, assay = "Spatial")

  image_cells <- c("cell3", "cell1", "cell4")
  coordinates <- data.frame(
    tissue = 1L,
    row = c(3L, 1L, 4L),
    col = c(13L, 11L, 14L),
    imagerow = c(30, 10, 40),
    imagecol = c(300, 100, 400),
    row.names = image_cells
  )
  scale_factors <- structure(
    list(spot = 4, fiducial = 8, hires = 0.25, lowres = 0.5),
    class = "scalefactors"
  )
  image <- methods::new(
    Class = "VisiumV1",
    assay = "Spatial",
    key = "slice1_",
    image = array(0.5, dim = c(6, 8, 3)),
    scale.factors = scale_factors,
    coordinates = coordinates
  )
  methods::slot(image, "spot.radius") <- Seurat::Radius(image, scale = "lowres")
  suppressWarnings(object[["slice1"]] <- image)

  embeddings <- matrix(
    seq_len(8),
    nrow = 4,
    dimnames = list(cells, c("pca_1", "pca_2"))
  )
  object[["pca"]] <- Seurat::CreateDimReducObject(
    embeddings = embeddings,
    key = "pca_",
    assay = "Spatial"
  )

  path <- tempfile(fileext = ".h5ad")
  on.exit(unlink(path), add = TRUE)

  expect_message(
    WriteH5AD_Spatial(
      seurat_object = object,
      output_path = path,
      library_id = "library1"
    ),
    "Spatial H5AD saved"
  )

  anndata <- reticulate::import("anndata", convert = FALSE)
  adata <- anndata$read_h5ad(path)
  obs_names <- as.character(reticulate::py_to_r(adata$obs_names$to_list()))
  spatial <- reticulate::py_to_r(adata$obsm[["spatial"]])
  pca <- reticulate::py_to_r(adata$obsm[["X_pca"]])
  spatial_uns <- reticulate::py_to_r(adata$uns[["spatial"]])

  # The image omits cell2, so every exported component must use the remaining
  # assay order rather than the original coordinate order.
  expect_equal(obs_names, c("cell1", "cell3", "cell4"))
  expect_equal(spatial[, 1], c(100, 300, 400))
  expect_equal(spatial[, 2], c(10, 30, 40))
  expect_equal(pca, unname(embeddings[c("cell1", "cell3", "cell4"), ]))

  library <- spatial_uns[["library1"]]
  expect_equal(dim(library$images$lowres), c(6L, 8L, 3L))
  expect_false("hires" %in% names(library$images))
  expect_equal(library$scalefactors$spot_diameter_fullres, 4)
  expect_equal(library$scalefactors$tissue_lowres_scalef, 0.5)

  roundtrip <- ReadH5AD(path)
  roundtrip_coordinates <- methods::slot(
    roundtrip[["slice1"]],
    "coordinates"
  )
  expect_equal(rownames(roundtrip_coordinates), obs_names)
  expect_equal(roundtrip_coordinates$imagerow, c(10, 30, 40))
  expect_equal(roundtrip_coordinates$imagecol, c(100, 300, 400))
})
