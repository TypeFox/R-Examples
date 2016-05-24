#' Wrapper for message and sprintf.
#'
#' A simple wrapper for \code{message(sprintf(...))}.
#'
#' @param ... [any]\cr
#'   See \code{\link{sprintf}}.
#' @param .newline [logical(1)]\cr
#'   Add a newline to the message. Default is \code{TRUE}.
#' @return Nothing.
#' @export
#' @examples
#' msg = "a message"
#' warningf("this is %s", msg)
messagef = function(..., .newline = TRUE) {
  message(sprintf(...), appendLF = .newline)
}
