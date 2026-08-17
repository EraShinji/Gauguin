test_that("test for reading anndata and parsing object rightly", {
  skip_if(!reticulate::py_module_available("anndata"), "anndata is not available")
  sce <- system.file("extdata", "sce.h5ad", package = "Gauguin")
  object <- ReadH5AD(sce)

  expect_s4_class(object, "Seurat")
  expect_gt(ncol(object), 0)
  expect_gt(nrow(object), 0)
})

test_that("ReadH5AD preserves Visium images, scale factors, and coordinate orientation", {
  skip_if(!reticulate::py_module_available("anndata"), "anndata is not available")

  path <- tempfile(fileext = ".h5ad")
  on.exit(unlink(path), add = TRUE)
  python_path <- normalizePath(path, winslash = "/", mustWork = FALSE)

  reticulate::py_run_string(paste0(
    "import anndata as ad\n",
    "import numpy as np\n",
    "import pandas as pd\n",
    "adata = ad.AnnData(\n",
    "    X=np.array([[1, 0], [0, 2], [3, 4]], dtype=np.int32),\n",
    "    obs=pd.DataFrame({\n",
    "        'in_tissue': [1, 1, 1],\n",
    "        'array_row': [0, 1, 2],\n",
    "        'array_col': [2, 3, 4]\n",
    "    }, index=['cell1', 'cell2', 'cell3']),\n",
    "    var=pd.DataFrame(index=['gene1', 'gene2'])\n",
    ")\n",
    "adata.obsm['spatial'] = np.array([[10, 20], [30, 40], [50, 60]], dtype=float)\n",
    "adata.obsm['c2l_cell_abundance_fraction'] = pd.DataFrame(\n",
    "    {'B cells': [0.1, 0.2, 0.3], 'T cells': [0.9, 0.8, 0.7]},\n",
    "    index=adata.obs_names\n",
    ")\n",
    "adata.uns['spatial'] = {\n",
    "    'library1': {\n",
    "        'images': {\n",
    "            'hires': np.zeros((12, 16, 3), dtype=float),\n",
    "            'lowres': np.zeros((6, 8, 3), dtype=float)\n",
    "        },\n",
    "        'scalefactors': {\n",
    "            'spot_diameter_fullres': 4.0,\n",
    "            'fiducial_diameter_fullres': 8.0,\n",
    "            'tissue_hires_scalef': 0.25,\n",
    "            'tissue_lowres_scalef': 0.5\n",
    "        }\n",
    "    }\n",
    "}\n",
    "adata.write_h5ad(r'", python_path, "')\n"
  ))

  object <- ReadH5AD(path)
  image <- object[["slice1"]]
  coordinates <- methods::slot(image, "coordinates")
  scale_factors <- Seurat::ScaleFactors(image)

  expect_equal(dim(methods::slot(image, "image")), c(6L, 8L, 3L))
  expect_equal(unname(scale_factors$spot), 4)
  expect_equal(unname(scale_factors$fiducial), 8)
  expect_equal(unname(scale_factors$hires), 0.25)
  expect_equal(unname(scale_factors$lowres), 0.5)
  expect_equal(coordinates$imagerow, c(20, 40, 60))
  expect_equal(coordinates$imagecol, c(10, 30, 50))
  expect_equal(methods::slot(image, "spot.radius"), 0.25)

  abundance <- Seurat::Embeddings(object, reduction = "c2l_cell_abundance_fraction")
  expect_true(is.matrix(abundance))
  expect_equal(dim(abundance), c(3L, 2L))
  expect_equal(
    colnames(abundance),
    c("c2lcellabundancefraction_1", "c2lcellabundancefraction_2")
  )
})
