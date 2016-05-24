#########################################################################/**
# @RdocDefault setCachePath
#
# @title "Sets the path to the file cache directory"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{dirs}{A @character @vector constituting the path to the
#      cache subdirectory of interest.}
#   \item{path}{The path to override the path according to the 
#      \code{dirs} argument.}
#   \item{...}{Not used.}
# }
#
# \value{
#   Returns nothing.
# }
#
# @author
#
# \seealso{
#   @see "getCachePath".
# }
#
# @keyword "programming"
# @keyword "IO"
# @keyword "internal"
#*/######################################################################### 
setMethodS3("setCachePath", "default", function(dirs=NULL, path=NULL, ...) {
  subname <- paste(dirs, collapse="/");
  name <- paste("R.cache:cachePath", subname, sep=":");
  opts <- list(path);
  names(opts) <- name;
  ovalue <- options(opts);
  invisible(ovalue);
}, export=FALSE)


############################################################################
# HISTORY:
# 2007-01-24
# o Created.
############################################################################
