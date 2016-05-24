#########################################################################/**
# @RdocDefault memoizedCall
#
# @title "Calls a function with memoization"
#
# \description{
#  @get "title", that is, caches the results to be retrieved if
#  the function is called again with the exact same arguments.
# }
#
# @synopsis
#
# \arguments{
#   \item{what}{The @function to be called, or a @character string
#     specifying the name of the function to be called,
#     cf. @see "base::do.call".}
#   \item{...}{Arguments passed to the function.}
#   \item{envir}{The @environment in which the function is evaluated.}
#   \item{force}{If @TRUE, any cached results are ignored, otherwise not.}
#   \item{sources, dirs}{Optional arguments passed to
#     @see "loadCache" and @see "saveCache".}
# }
#
# \value{
#   Returns the result of the function call.
# }
#
# \details{
#   If the @function returns @NULL, that particular function call is
#   \emph{not} memoized.
# }
#
# @author
#
# \seealso{
#  Internally, @see "loadCache" is used to load memoized results,
#  if available.  If not available, then @see "do.call" is used to
#  evaluate the function call,
#  and @see "saveCache" is used to save the results to cache.
# }
#
# @keyword "programming"
# @keyword "IO"
#*/#########################################################################
setMethodS3("memoizedCall", "default", function(what, ..., envir=parent.frame(), force=FALSE, sources=NULL, dirs=NULL) {
  # 1. Generate cache file
  key <- list(what=what, ...);
  pathnameC <- generateCache(key=key, dirs=dirs);

  # 1. Look for memoized results
  if (!force) {
    res <- loadCache(pathname=pathnameC, sources=sources);
    if (!is.null(res)) return(res)
  }

  # 2. Otherwise, call method with arguments
  res <- do.call(what, args=list(...), quote=FALSE, envir=envir);

  # 3. Memoize results
  saveCache(res, pathname=pathnameC, sources=sources);

  # 4. Return results
  res;
}) # memoizedCall()



#######################################################################
# HISTORY:
# 2015-02-27
# o SPEEDUP: Now memoizedCall() generates cache pathname only once.
# 2011-02-14
# o Added argument 'sources' to memoizedCall().
# 2011-02-13
# o Created.
#######################################################################
