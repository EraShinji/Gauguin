#' ReadH5AD
#' @description
#' Read H5AD AnnData Object,and convert to Seurat Object
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
#' @importFrom methods as
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

  # 判断标准：是否存在明显的负值
  is_scaled = any(as.numeric(sample_X) < -0.0001)
  if ("scaled" %in% layers_list) {
    scaled_mtx = t(py_to_r(adata$layers["scaled"]))
    rownames(scaled_mtx) = as.character(py_to_r(adata$var_names$to_list()))
    colnames(scaled_mtx) = obs_names
    seurat_obj = SetAssayData(seurat_obj, slot = "scale.data", new.data = as.matrix(scaled_mtx))
  } else if(is_scaled){
    message("adata.X might be scaled, please consider use adata.X as layer `scale.data`")
  }else{
    message("No scaled.data be yieled in anndata object")
  }

  obsm_dict = py_call(adata$obsm$as_dict)

  for (key in names(obsm_dict)) {
    embedding = obsm_dict[[key]]
    embedding= py_to_r(embedding)
    rownames(embedding) = colnames(seurat_obj)
    colnames(embedding) = paste0(key, "_", 1:ncol(embedding))
    fixed_key = gsub("^X_", "", key)
    fixed_key = paste0(fixed_key, "")
    seurat_obj[[fixed_key]] = CreateDimReducObject(embeddings = embedding,
                                                    key = fixed_key, assay = "RNA")
  }
  return(seurat_obj)
  }
