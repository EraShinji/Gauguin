test_that("test for saving anndata and parsing object rightly", {
  sce <- readRDS(system.file("extdata", "sce.rds", package = "Gauguin"))
  WriteH5AD(sce,"C:/Users/aleclanned/.local/share/mamba/envs/Single.Cell.Analysis","sce2.h5ad")
  }
)
test_that("test for saving anndata SCTransform and parsing object rightly", {
  sce <- readRDS(system.file("extdata", "sce_sct.rds", package = "Gauguin"))
  WriteH5AD(sce,"C:/Users/aleclanned/.local/share/mamba/envs/Single.Cell.Analysis","sce2.h5ad",assay = "SCT")
}
)
