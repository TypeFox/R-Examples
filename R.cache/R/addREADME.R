.addREADME <- function(to=getCacheRootPath(), ...) {
  # Add a README.txt to cache root (expaining what the directory is)
  filename <- "README.txt";
  pathnameD <- file.path(to, filename);
  if (!isFile(pathnameD)) {
    pathnameS <- system.file("_Rcache", filename, package="R.cache");
    file.copy(pathnameS, pathnameD);
  }
} # .addREADME()
 

############################################################################
# HISTORY:
# 2012-11-28
# o Added internal .addREADME().
# o Created.
############################################################################
