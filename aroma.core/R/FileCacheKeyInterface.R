###########################################################################/**
# @RdocClass FileCacheKeyInterface
#
# @title "The FileCacheKeyInterface class interface"
#
# \description{
#  @classhierarchy
# }
# 
# @synopsis
#
# \arguments{
#   \item{...}{Not used.}
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"
# }
#
# @author
#
# @keyword internal
#*/########################################################################### 
setConstructorS3("FileCacheKeyInterface", function(...) {
  extend(CacheKeyInterface(), "FileCacheKeyInterface");
})



###########################################################################/**
# @RdocMethod getCacheKey
#
# @title "Gets a list of cache key items"
#
# \description{
#  @get "title" that will be added to other cache key items used to 
#  generate the cache key.
# }
#
# @synopsis
#
# \arguments{
#  \item{...}{Additional named arguments to be appended to the list
#    of key items.}
# }
#
# \value{
#  Returns a @list of cache items.
# }
#
# \details{
#  The default list of cache key items are:
#  \itemize{
#   \item the class name of the object as a @character string.
#   \item the file name (aka the basename) as a @character string.
#   \item the file size as a @numeric scalar.
#   \item the file checksum as a @character string.
#  }
#
#  Classes extending/implementing this @see "R.oo::Interface" may override 
#  these items.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/###########################################################################
setMethodS3("getCacheKey", "FileCacheKeyInterface", function(this, ...) {
  keys <- NextMethod("getCacheKey");
  keys$filename <- getFilename(this);
  keys$filesize <- getFileSize(this);
  keys$checksum <- getChecksum(this);
  keys;
})



############################################################################
# HISTORY:
# 2012-11-28
# o Added FileCacheKeyInterface
# o Created.
############################################################################
