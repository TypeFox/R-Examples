#########################################################################/**
# @RdocDefault generateCache
#
# @title "Generates a cache pathname from a key object"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{key}{A @list or an @environment from which a @character string
#     checksum will be calculated and that will constitute the name part
#     of the cache filename.}
#   \item{suffix}{A @character string to be appended to the end of the
#     filename.}
#   \item{...}{Arguments passed to @see "getCachePath".}
# }
#
# \value{
#   Returns the pathname as a @character string.
# }
#
# @author
#
# \seealso{
#   Internally, the generic function @see "getChecksum" is used to
#   calculate the checksum of argument \code{key}.
# }
#
# @keyword "programming"
# @keyword "IO"
# @keyword "internal"
#*/#########################################################################
setMethodS3("generateCache", "default", function(key, suffix=".Rcache", ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'key':
  if (!is.list(key) && !is.environment(key)) {
    throw("Argument 'key' must be a list, an environment: ", class(key)[1L]);
  }

  # Generate cache name from hash code of key object
  cacheName <- getChecksum(key);

  # Add cache directory or pathname
  path <- getCachePath(...);
  cacheName <- file.path(path, cacheName);

  # Add suffix
  cacheName <- paste(cacheName, suffix, sep="");

  cacheName;
}, export=FALSE)


############################################################################
# HISTORY:
# 2014-01-13
# o BUG FIX: generateCache(key=NULL) would generate "Error in
#   file.path(path, cacheName) : object 'cacheName' not found".
# 2011-04-01
# o Now generateCache() utilizes a generic function getChecksum()
#   to obtain the checksum.
# 2007-07-02
# o Added support for argument 'key' is an environment.
# 2005-12-09
# o Removed 'file' argument.
# 2005-12-06
# o Created.
############################################################################
