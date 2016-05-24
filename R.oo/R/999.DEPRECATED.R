## covr: skip=all
###########################################################################/**
# @set class=Object
# @RdocMethod gc
#
# @title "Clear cached fields and calls the garbage collector"
#
# \description{
#  @get "title".  Cached fields are set to @NULL when cleared.
#
#  \emph{This method is deprecated.
#   Please use \code{clearCache(..., gc=TRUE)} instead.}
# }
#
# @synopsis
#
# \arguments{
#   \item{...}{Not used.}
# }
#
# \value{
#   Returns nothing.
# }
#
# @examples "../incl/gc.Object.Rex"
#
# @author
#
# \seealso{
#   To clear the fields without calling the garbage collector,
#   see @seemethod "clearCache".
#   @seeclass
# }
#
# @keyword programming
# @keyword methods
#*/###########################################################################
setMethodS3("gc", "Object", function(this, ...) {
  .Deprecated(msg="Use clearCache(..., gc=TRUE) instead.");
  clearCache(this, gc=TRUE);
}, deprecated=TRUE)


###########################################################################/**
# @set class=Object
# @RdocMethod registerFinalizer
#
# @title "Registers a finalizer hook for the object [DEFUNCT]"
#
# \description{
#  @get "title".
#  The finalizer hook calls @seemethod "finalize" on the @see Object when
#  it is garbage collected.
#  This method is only intended to be called inside the constructor, if
#  at all.
# }
#
# @synopsis
#
# \arguments{
#   \item{...}{Not used.}
# }
#
# \value{
#  Returns nothing.
# }
#
# @author
#
# \seealso{
#   Internally, @see "base::reg.finalizer" is used.
#   @seeclass
# }
#
# @keyword programming
# @keyword methods
# @keyword internal
#*/###########################################################################
setMethodS3("registerFinalizer", "Object", function(this, ...) {
  .Defunct("registerFinalizer() for Object is deprecated.");
}, protected=TRUE, deprecated=TRUE) # registerFinalizer()


#########################################################################/**
# @set class=Package
# @RdocMethod update
#
# @title "Updates the package is a newer version is available"
#
# \description{
#   \emph{This method is defunct. Use @see "utils::update.packages" instead.}
# }
#
# @synopsis
#
# \arguments{
#   \item{contribUrl}{The URL from where the package can be installed and
#    updated. By default the URL according to the DESCRIPTION is assumed.
#    If the URL is missing, CRAN is assumed.}
#   \item{force}{If @TRUE, the package will reinstalled even if it is
#    up to date according to the version number.}
#   \item{verbose}{If @TRUE, more detailed information is returned.}
#   \item{...}{Not used.}
# }
#
# \value{
#   Returns (invisibly) @TRUE if the package was updated, otherwise @FALSE.
# }
#
# @author
#
# \seealso{
#   @see "utils::update.packages".
#   @seeclass
# }
#
# @keyword internal
#*/#########################################################################
setMethodS3("update", "Package", function(object, contribUrl=c(getContribUrl(object), getDevelUrl(object)), force=FALSE, reload=TRUE, verbose=TRUE, ...) {
  .Defunct(msg=sprintf("update() for Package is defunct. Use update.packages(\"%s\") instead."), getName(object));
}, protected=TRUE, deprecated=TRUE)


############################################################################
# HISTORY:
# 2015-01-05
# o CLEANUP: Defunct update() for Package.
# 2013-09-25
# o CLEANUP: Deprecated update() for Package.
# 2014-02-22
# o DEPRECATED: Deprecated gc() for Object.  Use clearCache(..., gc=TRUE)
#   instead.
# 2014-01-05
# o CLEANUP: Defunct registerFinalizer() for Object.
# o Created.
############################################################################
