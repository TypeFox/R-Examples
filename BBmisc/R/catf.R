#' Wrapper for cat and sprintf.
#'
#' A simple wrapper for \code{cat(sprintf(...))}.
#'
#' @param ... [any]\cr
#'   See \code{\link{sprintf}}.
#' @param file [\code{character(1)}]\cr
#'   See \code{\link{cat}}.
#'   Default is \dQuote{}.
#' @param append [\code{logical(1)}]\cr
#'   See \code{\link{cat}}.
#'   Default is \code{FALSE}.
#' @param newline [\code{logical(1)}]\cr
#'   Append newline at the end?
#'   Default is \code{TRUE}.
#' @return Nothing.
#' @export
#' @examples
#' msg = "a message."
#' catf("This is %s", msg)
catf = function(..., file = "", append = FALSE, newline = TRUE) {
  cat(sprintf(...), ifelse(newline, "\n", ""), sep = "", file = file, append = append)
}
