#' Wrapper for warning and sprintf.
#'
#' A wrapper for \code{\link{warning}} with \code{\link{sprintf}} applied to the arguments.
#'
#' @param ... [any]\cr
#'   See \code{\link{sprintf}}.
#' @param immediate [\code{logical(1)}]\cr
#'   See \code{\link{warning}}.
#'   Default is \code{TRUE}.
#' @param warning.length [\code{integer(1)}]\cr
#'   Number of chars after which the warning message
#'   gets truncated, see ?options.
#'   Default is 8170.
#' @return Nothing.
#' @export
#' @examples
#' msg = "a warning"
#' warningf("this is %s", msg)
warningf = function(..., immediate = TRUE, warning.length = 8170L) {
  msg = sprintf(...)
  if (immediate) {
    old = getOption("warn")
    # dont change warn setting if it is 2 (= error)
    if (old <= 0L) {
      on.exit(options(warn = old))
      options(warn = 1L)
    }
  }
  obj = simpleWarning(msg, call = sys.call(sys.parent()))
  warning(obj)
}
