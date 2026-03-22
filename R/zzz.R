.onAttach <- function(libname, pkgname) {
  version <- utils::packageVersion(pkgname)
  msg <- paste0(
    "Gauguin (v", version, ")\n",
    "ℹ️  To get started, visit the tutorials: https://www.geneinnova.dev/gauguin\n",
    "   Or run: browseVignettes('Gauguin')"
  )
  packageStartupMessage(msg)
}
