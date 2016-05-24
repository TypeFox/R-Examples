#' Extract wrapped function.
#'
#' The \pkg{smoof} package offers means to let a function log its evaluations or
#' even store to points and function values it has been evaluated on. This is done
#' by wrapping the function with other functions. This helper function extract
#' the wrapped function.
#'
#' @note If this function is applied to a simple \code{smoof_function}, the
#' \code{smoof_function} itself is returned.
#'
#' @param fn [\code{smoof_wrapped_function}]\cr
#'   Wrapping function.
#' @param deepest [\code{logical(1)}]\cr
#'   Function may be wrapped with multiple wrappers. If \code{deepest} is set to
#'   \code{TRUE} the function unwraps recursively until the \dQuote{deepest} wrapped
#'   \code{smoof_function} is reached. Default is \code{FALSE}.
#' @return [\code{function}]
#' @seealso \code{\link{addCountingWrapper}}, \code{\link{addLoggingWrapper}}
#' @export
getWrappedFunction = function(fn, deepest = FALSE) {
  UseMethod("getWrappedFunction")
}

#' @export
getWrappedFunction.smoof_wrapped_function = function(fn, deepest = FALSE) {
  wrapped.fn = environment(fn)$fn
  if (deepest) {
    return(getWrappedFunction(wrapped.fn, deepest))
  }
  return(wrapped.fn)
}

#' @export
getWrappedFunction.smoof_function = function(fn, deepest = FALSE) {
  return(fn)
}
