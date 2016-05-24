#' A caching wrapper around load2.
#'
#' This closure returns a wrapper around \code{\link{load2}} which per
#' default caches loaded objects and returns the cached version
#' in subsequent calls.
#'
#' @param use.cache [\code{logical(1)}]\cr
#'  Enable the cache?
#'  Default is \code{TRUE}.
#' @return [\code{function()}] with argument \code{slot}
#'  (name of the slot to cache the object in, default is \dQuote{default}).
#'  All other arguments are passed down to \code{\link{load2}}.
#' @export
makeFileCache = function(use.cache = TRUE) {
  assertFlag(use.cache)
  .cache = list()

  function(file, slot = "default", ...) {
    if (use.cache) {
      if (is.null(.cache[[slot]]) || .cache[[slot]]$file != file)
        .cache[[slot]] = list(file = file, obj = load2(file = file, ...))
      return(.cache[[slot]]$obj)
    }
    return(load2(file = file, ...))
  }
}
