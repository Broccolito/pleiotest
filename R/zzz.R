.onLoad = function(libname, pkgname){
  if (!requireNamespace("BiocManager", quietly = TRUE)){
    message("Installing metaCCA from Bioconductor...")
    BiocManager::install("metaCCA", ask = FALSE)
  }
  if (!requireNamespace("metaCCA", quietly = TRUE)){
    message("Installing metaCCA from Bioconductor...")
    BiocManager::install("metaCCA", ask = FALSE)
  }
  if (!requireNamespace("pleiosim", quietly = TRUE)){
    message("Installing pleiosim from GitHub")
    devtools::install_github("Broccolito/pleiosim")
  }
  if (!requireNamespace("MPAT", quietly = TRUE)){
    message("Installing MPAT from GitHub")
    devtools::install_github("Broccolito/MPATclone")
  }
  if (!requireNamespace("MTARclone", quietly = TRUE)){
    message("Installing MTAR from GitHub")
    devtools::install_github("Broccolito/MTARclone")
  }

}

