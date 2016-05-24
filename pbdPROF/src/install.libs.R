### Modified from Rserve/src/install.libs.R
### For libs
files <- c("pbdPROF.so", "pbdPROF.so.dSYM", "pbdPROF.dylib", "pbdPROF.dll",
           "symbols.rds")
files <- files[file.exists(files)]
if(length(files) > 0){
  libsarch <- if (nzchar(R_ARCH)) paste("libs", R_ARCH, sep='') else "libs"
  dest <- file.path(R_PACKAGE_DIR, libsarch)
  dir.create(dest, recursive = TRUE, showWarnings = FALSE)
  file.copy(files, dest, overwrite = TRUE, recursive = TRUE)
}
### For lib
file <- "libfpmpi.a"
if(file.exists(file)){  ### External libraries do not build "libfpmpi.a"
  libarch <- if (nzchar(R_ARCH)) paste("lib", R_ARCH, sep='') else "lib"
  dest <- file.path(R_PACKAGE_DIR, libarch)
  dir.create(dest, recursive = TRUE, showWarnings = FALSE)
  file.copy(file, dest, overwrite = TRUE)
}
### For etc
file <- "Makeconf"
if(file.exists(file)){
  etcarch <- if (nzchar(R_ARCH)) paste("etc", R_ARCH, sep='') else "etc"
  dest <- file.path(R_PACKAGE_DIR, etcarch)
  dir.create(dest, recursive = TRUE, showWarnings = FALSE)
  file.copy(file, dest, overwrite = TRUE)
}
