# Gauguin
**G**eneralized **A**nndata **U**tility for **G**uaranteed **U**ninterrupted **I**ntegration into **N**ative-Seurat
 Generalized Anndata Utility for Guaranteed Uninterrupted Integration into Native-Seurat

Gauguin is a high-performance bridge for single-cell data analysis. It enables the seamless transfer of genomic landscapes from Python's AnnData (H5AD) objects to R's Seurat ecosystem, ensuring your workflow remains uninterrupted across different programming environments.

## Key Features
- Zero-Friction Integration: Effortlessly convert .h5ad files into native Seurat objects.

- Minimal Dependency: Designed to be lightweight. You don't need a bloated environment; as long as Seurat is present, Gauguin handles the rest.

- Full Data Preservation: Automatically migrates raw counts, normalized layers, and dimensional embeddings (e.g., PCA, UMAP, t-SNE) while maintaining metadata integrity.

- Native Performance: Built to handle large-scale single-cell datasets with high efficiency。

## Install
You can install the development version of Gauguin directly from GitHub:
```
# Install devtools if you haven't already
if (!require("devtools", quietly = TRUE))
    install.packages("devtools")

# Install Gauguin
devtools::install_github("EraShinji/Gauguin")
```
Alternatively, you can install the package locally:
```
install.packages("path/to/Gauguin", repos = NULL, type = "source")
```
To get started, visit the tutorials: https://www.geneinnova.dev/gauguin
We welcome feedback and contributions! If you encounter any bugs or have feature requests, please open an issue on GitHub.
