#' Wrapper for stop and sprintf.
#'
#' A wrapper for \code{\link{stop}} with \code{\link{sprintf}} applied to the arguments.
#' Notable difference is that error messages are not truncated to 1000 characters
#' by default.
#'
#' @param ... [any]\cr
#'   See \code{\link{sprintf}}.
#' @param warning.length [\code{integer(1)}]\cr
#'   Number of chars after which the error message
#'   gets truncated, see ?options.
#'   Default is 8170.
#' @return Nothing.
#' @export
#' @examples
#' err = "an error."
#' try(stopf("This is %s", err))
stopf = function(..., warning.length = 8170L) {
  msg = sprintf(...)
  obj = simpleError(msg, call = sys.call(sys.parent()))
  old.opt = getOption("warning.length")
  on.exit(options(warning.length = old.opt))
  options(warning.length = warning.length)
  stop(obj)
}
