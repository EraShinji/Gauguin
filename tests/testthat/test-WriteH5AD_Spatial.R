# test-WriteH5AD_Spatial.R
# Uses the stxBrain (10x Visium anterior1) dataset from SeuratData,
# following the official Seurat spatial vignette.

# --- Setup: load stxBrain anterior1 once for all tests ---
skip_if_not_installed("SeuratData")

library(SeuratData)
brain <- tryCatch(
  LoadData("stxBrain", type = "anterior1"),
  error = function(e) NULL
)

# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------

test_that("WriteH5AD_Spatial writes a valid h5ad file for stxBrain VisiumV1", {
  skip_on_cran()
  skip_if(is.null(brain), "stxBrain dataset not available")
  out_path <- tempfile(fileext = ".h5ad")
  on.exit(unlink(out_path), add = TRUE)

  expect_no_error(
    WriteH5AD_Spatial(
      seurat_object   = brain,
      output_path     = out_path,
      assay           = "Spatial",
      library_id      = "anterior1",
      image_resolution = "lowres"
    )
  )
  expect_true(file.exists(out_path))
  expect_gt(file.size(out_path), 0)
})

test_that("WriteH5AD_Spatial preserves lowres when both resolutions are requested", {
  skip_on_cran()
  skip_if(is.null(brain), "stxBrain dataset not available")
  out_path <- tempfile(fileext = ".h5ad")
  on.exit(unlink(out_path), add = TRUE)

  expect_warning(
    WriteH5AD_Spatial(
      seurat_object    = brain,
      output_path      = out_path,
      image_resolution = "both"
    ),
    "no hires image is available"
  )
  expect_true(file.exists(out_path))
})

test_that("WriteH5AD_Spatial falls back when specified assay is missing", {
  skip_on_cran()
  skip_if(is.null(brain), "stxBrain dataset not available")
  out_path <- tempfile(fileext = ".h5ad")
  on.exit(unlink(out_path), add = TRUE)

  expect_message(
    WriteH5AD_Spatial(
      seurat_object = brain,
      output_path   = out_path,
      assay         = "NonExistentAssay"
    ),
    "not found|Using"
  )
  expect_true(file.exists(out_path))
})

test_that("WriteH5AD_Spatial errors on invalid image_name", {
  skip_on_cran()
  skip_if(is.null(brain), "stxBrain dataset not available")
  out_path <- tempfile(fileext = ".h5ad")
  on.exit(unlink(out_path), add = TRUE)

  expect_error(
    WriteH5AD_Spatial(
      seurat_object = brain,
      output_path   = out_path,
      image_name    = "nonexistent_slice"
    ),
    "not found"
  )
})

test_that("WriteH5AD_Spatial warns and writes non-spatial h5ad for plain object", {
  skip_on_cran()
  # Build a minimal non-spatial Seurat object
  counts <- Matrix::Matrix(
    matrix(rpois(50, 5), nrow = 5, ncol = 10,
           dimnames = list(paste0("G", 1:5), paste0("C", 1:10))),
    sparse = TRUE
  )
  plain_obj <- Seurat::CreateSeuratObject(counts = counts, assay = "Spatial")
  out_path  <- tempfile(fileext = ".h5ad")
  on.exit(unlink(out_path), add = TRUE)

  expect_warning(
    WriteH5AD_Spatial(
      seurat_object = plain_obj,
      output_path   = out_path
    ),
    "No spatial images"
  )
  expect_true(file.exists(out_path))
})

test_that("WriteH5AD_Spatial output roundtrips correctly via anndata", {
  skip_on_cran()
  skip_if(is.null(brain), "stxBrain dataset not available")
  out_path <- tempfile(fileext = ".h5ad")
  on.exit(unlink(out_path), add = TRUE)

  WriteH5AD_Spatial(
    seurat_object = brain,
    output_path   = out_path,
    library_id    = "anterior1"
  )

  anndata <- reticulate::import("anndata", delay_load = FALSE, convert = FALSE)
  adata   <- anndata$read_h5ad(out_path)

  n_cells <- ncol(Seurat::GetAssayData(brain, assay = "Spatial", layer = "counts"))
  n_genes <- nrow(Seurat::GetAssayData(brain, assay = "Spatial", layer = "counts"))

  # Dimensions

  expect_equal(as.integer(reticulate::py_to_r(adata$n_obs)), n_cells)
  expect_equal(as.integer(reticulate::py_to_r(adata$n_vars)), n_genes)

  # Spatial coordinates written to obsm
  builtins <- reticulate::import_builtins(convert = FALSE)
  obsm_keys <- as.character(reticulate::py_to_r(builtins$list(adata$obsm$keys())))
  expect_true("spatial" %in% obsm_keys)
  spatial_mat <- reticulate::py_to_r(adata$obsm[["spatial"]])
  expect_equal(nrow(spatial_mat), n_cells)
  expect_equal(ncol(spatial_mat), 2L)

  # uns["spatial"] structure
  uns <- reticulate::py_to_r(adata$uns)
  expect_true("spatial" %in% names(uns))
  spatial_uns <- uns[["spatial"]]
  expect_true("anterior1" %in% names(spatial_uns))
  expect_true("images"       %in% names(spatial_uns[["anterior1"]]))
  expect_true("scalefactors" %in% names(spatial_uns[["anterior1"]]))
})
