# This script was modified from the pdbBASE package.
#
files <- c("libboom.a", "symbols.rds")
files <- files[file.exists(files)]
if (length(files) > 0) {
  if (nzchar(R_ARCH)) {
      libsarch <- paste("lib", R_ARCH, sep='')
  } else {
    libsarch <- "lib"
  }
  dest <- file.path(R_PACKAGE_DIR, libsarch)
  dir.create(dest, recursive = TRUE, showWarnings = FALSE)
  file.copy(files, dest, overwrite = TRUE, recursive = TRUE)
}
