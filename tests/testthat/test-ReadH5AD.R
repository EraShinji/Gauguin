test_that("test for reading anndata and parsing object rightly", {
  sce <- system.file("extdata", "sce.h5ad", package = "Gauguin")
  ReadH5AD(sce,"C:/Users/aleclanned/.local/share/mamba/envs/Single.Cell.Analysis")
})
