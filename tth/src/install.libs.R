binaries <- c("tth", "ttm")
if(WINDOWS) binaries <- paste(binaries, ".exe", sep="")

if(all(file.exists(binaries))) {
  libarch <- if (nzchar(R_ARCH)) paste('libs', R_ARCH, sep='') else 'libs'
  dest <- file.path(R_PACKAGE_DIR, libarch)
  dir.create(dest, recursive = TRUE, showWarnings = FALSE)
  file.copy(binaries, dest, overwrite = TRUE)
}
