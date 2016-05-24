.setupCacheRootPath <- function(...) {
  # Setup the cache root path, possibly by prompting the user.
  ns <- getNamespace("R.cache");
  setupCacheRootPath <- get("setupCacheRootPath", mode="function", envir=ns);
  setupCacheRootPath();
} # .setupCacheRootPath()

# CRAN POLICY: Add precalculated memoization files to the R.cache
# directory, unless running interactively.  The reason for doing this
# is solely to make segmentBy[Non]PairedPSCBS examples to run faster
# on R CMD check but not having to create these memoized files.
# /HB 2012-11-05
# UPDATE: Now it will also gain first-time users. /HB 2013-09-27
.prememoize <- function(verbose=FALSE) {
  # Explictly setup cache root here, since it's only done by 'R.cache'
  # if that package is attached.  Here we only load it. /HB 2013-09-27
  .setupCacheRootPath();

  # This will make sure that the pre-generated calculations available
  # in the 'PSCBS' package are copied to the R.cache cache directory.
  # This regardless of whether a 'PSCBS' cache subdirectory exists
  # or not. /HB 2013-09-27
  path <- "PSCBS/segmentByCBS/sbdry"
  pathS <- system.file("misc/_Rcache", path, package="PSCBS");
  pathD <- getCachePath(path);
  copyDirectory(pathS, pathD, copy.mode=FALSE, recursive=FALSE, overwrite=TRUE);
  if (verbose) {
    message("Added pre-memoized calculations: ", getAbsolutePath(pathD));
  }
} # .prememoize()

############################################################################
# HISTORY:
# 2013-09-27
# o Now .prememorize() also copies pre-generated calculations in
#   interactive session.  It is also called every time the package
#   is attached, which means it will also gain first-time users.
# o Added .setupCacheRootPath() until R.cache exports it.
# o Added argument 'verbose' to .prememoize().
# 2013-09-26
# o CLEANUP: Now .prememoize() no longer attaches 'R.cache', but only
#   loads its namespace.
# 2012-11-05
# o Created.
############################################################################
