mydebug <- FALSE
if (mydebug) {
  cat("R_ARCH=", R_ARCH,"\n")
  cat("R_PACKAGE_DIR=", R_PACKAGE_DIR, "\n")
  cat("R_PACKAGE_NAME=", R_PACKAGE_NAME, "\n")
  cat("R_PACKAGE_SOURCE=", R_PACKAGE_SOURCE, "\n")
  cat("SHLIB_EXT=", SHLIB_EXT, "\n")
  cat("WINDOWS=", WINDOWS, "\n")
}
binary <- if(WINDOWS) "BayesX.exe" else "BayesX"
if ( file.exists(binary) ) {
  libarch <- if (nzchar(R_ARCH)) paste('libs', R_ARCH, sep='') else 'libs'
  dest <- file.path(R_PACKAGE_DIR, libarch)
  dir.create(dest, recursive = TRUE, showWarnings = FALSE)
  file.copy(binary, dest, overwrite = TRUE)
}

