#' WriteH5AD
#' @description
#' Write Seurat object into h5ad AnnData object format
#' @details
#' Warn: To keep the integrity and reproducibility, the `WriteH5AD` func will use counts layer as adata.X.
#'
#'
#' @param seurat_object The Seurat object required to be converted
#' @param env_path A path used to execute python environment.
#' @param output_path The path stored converted AnnData object
#' @param assay The Name of the assay to use for expression data. Default "RNA".
#'
#' @return No data needed to be return
#' @export

WriteH5AD = function(seurat_object,env_path,output_path,assay = "RNA"){

  if (!missing(env_path)) {
    use_python(env_path, required = TRUE)
  }

  tryCatch({
    anndata = reticulate::import("anndata", delay_load = FALSE)
    np = reticulate::import("numpy", delay_load = FALSE)
    pd = reticulate::import("pandas", delay_load = FALSE)
  }, error = function(e) {
    stop("Error: The required packages numpy, pandas or scanpy could not be found in the provided environment. Please check the environment.", e$message)
  })

  counts_matrix = Seurat::GetAssayData(seurat_obj, assay = assay, layer = 'counts')
  if (!inherits(counts_matrix, "matrix") && !inherits(counts_matrix, "dgCMatrix")) {
    counts_matrix = as.matrix(counts_matrix)
  }

  original_cell_ids = colnames(counts_matrix)
  gene_names = rownames(counts_matrix)
  counts_matrix = Matrix::t(counts_matrix)
  counts_matrix = as(counts_matrix, "dgCMatrix")


  meta_data = seurat_obj@meta.data
  meta_data$barcode = rownames(meta_data)

  adata = anndata$AnnData(X = counts_matrix, obs = pd$DataFrame(meta_data))
  adata$var_names = np$array(gene_names)

  if (length(seurat_obj@reductions) > 0) {
    for (reduction_name in names(seurat_obj@reductions)) {
      embeddings = Embeddings(seurat_obj, reduction = reduction_name)
      obsm_key = paste0("X_", reduction_name)
      adata$obsm[obsm_key] = np$array(embeddings)
    }
  } else {
    message("The embedding dimmed doesn't seem to be generated, pass.")
  }
  adata$write_h5ad(output_path)
  message(paste("Saved finished,","stored at",output_path))

}
