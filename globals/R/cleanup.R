#' @export
cleanup <- function(...) UseMethod("cleanup")

#' Drop certain types of globals
#'
#' @param globals A Globals object.
#' @param drop A character vector specifying what type of globals to drop.
#' @param \dots Not used
#'
#' @aliases cleanup
#' @export
cleanup.Globals <- function(globals, drop=c("base-packages"), ...) {
  where <- attr(globals, "where")

  names <- names(globals)
  keep <- rep(TRUE, times=length(globals))
  names(keep) <- names

  ## Drop objects that are part of one of the "base" packages
  if ("base-packages" %in% drop) {
    for (name in names) {
      if (isBasePkgs(environmentName(where[[name]]))) keep[name] <- FALSE
    }
  }

  ## Drop objects that are primitive functions
  if ("primitives" %in% drop) {
    for (name in names) {
      if (is.primitive(globals[[name]])) keep[name] <- FALSE
    }
  }

  ## Drop objects that calls .Internal()
  if ("internals" %in% drop) {
    for (name in names) {
      if (is.internal(globals[[name]])) keep[name] <- FALSE
    }
  }

  if (!all(keep)) {
    globals <- globals[keep]
  }

  globals
} # cleanup()
