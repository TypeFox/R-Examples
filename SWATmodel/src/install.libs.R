## Modified from Rserve
libarch = if (nzchar(version$arch)) paste("libs", version$arch, sep = "/") else "libs"
dest <- file.path(R_PACKAGE_DIR, libarch)
swats<- c("rswat2005.exe","rswat2009.exe","rswat2012.exe")
for(filename in swats) {
  dir.create(dest, recursive = TRUE, showWarnings = FALSE)
  file.copy(filename, dest, overwrite = TRUE)
} 

