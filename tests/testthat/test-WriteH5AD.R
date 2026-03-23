test_that("test for saving anndata and parsing object rightly", {
  test_sce <- readRDS("./inst/extdata/sce.rds")
  WriteH5AD(test_sce,"inst/extdata/sce2.h5ad","C:/Users/aleclanned/.local/share/mamba/envs/Single.Cell.Analysis")
  }
)
