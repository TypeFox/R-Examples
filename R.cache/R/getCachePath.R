#########################################################################/**
# @RdocDefault getCachePath
#
# @title "Gets the path to the file cache directory"
#
# \description{
#  @get "title".
#  If missing, the directory is created.
# }
#
# @synopsis
#
# \arguments{
#   \item{dirs}{A @character @vector constituting the path to the
#      cache subdirectory (of the \emph{cache root directory}
#      as returned by @see "getCacheRootPath") to be used.
#      If @NULL, the path will be the cache root path.}
#   \item{path, rootPath}{(Advanced) @character strings specifying the
#      explicit/default cache path and root cache path.}
#   \item{...}{Not used.}
# }
#
# \value{
#   Returns the path as a @character string.
# }
#
# @author
#
# \seealso{
#   @see "setCachePath".
# }
#
# @keyword "programming"
# @keyword "IO"
# @keyword "internal"
#*/#########################################################################
setMethodS3("getCachePath", "default", function(dirs=NULL, path=NULL, rootPath=getCacheRootPath(), ...) {
  # Get path where to store cache file
  if (is.null(path)) {
    # (1) Get/make default path
    # (a) Get path from options
    subname <- paste(dirs, collapse="/");
    name <- paste("R.cache:cachePath", subname, sep=":");
    path <- getOption(name);

    # (b) If not availble, make on
    path <- paste(c(rootPath, dirs), collapse=.Platform$file.sep);
  } else if (!isAbsolutePath(path)) {
    # (2) Get/make default path
    path <- file.path(rootPath, path);
  }

  # Create missing directory?
  if (!isDirectory(path)) {
    mkdirs(path);
    if (!isDirectory(path)) {
      throw("Could not create cache directory: ", path);
    }

    # Add a README.txt file, if missing.
    .addREADME(to=rootPath);
  }

  path;
}) # getCachePath()




############################################################################
# HISTORY:
# 2013-12-21
# o Added argument 'path' and 'rootPath' to getCachePath().
# 2012-09-10
# o Renamed the installed .Rcache/ directory to _Rcache/ to avoid
#   R CMD check NOTEs.
# 2011-12-29
# o Now getCachePath() adds a README.txt file to the root path, iff
#   missing. It explains why the directory structure exists and what
#   created it.
# 2007-01-24
# o Created.
############################################################################
