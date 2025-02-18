.onLoad <- function(libname, pkgname) {
  if (!requireNamespace("metaCCA", quietly = TRUE)) {
    message("Installing metaCCA from Bioconductor...")
    BiocManager::install("metaCCA", ask = FALSE)
    message("Installing pleiosim from GitHub")
    devtools::install_github("Broccolito/pleiosim")
    message("Installing MPAT from GitHub")
    devtools::install_github("Broccolito/MPATclone")
    message("Installing MTAR from GitHub")
    devtools::install_github("Broccolito/MTARclone")
  }
}

