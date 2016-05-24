#' Convert Integers to Strings
#'
#' This is the counterpart of \code{\link[base]{strtoi}}.
#' For a base greater than \sQuote{10}, letters \sQuote{a} to \sQuote{z}
#' are used to represent \sQuote{10} to \sQuote{35}.
#'
#' @param x [\code{integer}]\cr
#'  Vector of integers to convert.
#' @param base [\code{integer(1)}]\cr
#'  Base for conversion. Values between 2 and 36 (inclusive) are allowed.
#' @return \code{character(length(x))}.
#' @export
#' @examples
#' # binary representation of the first 10 natural numbers
#' itostr(1:10, 2)
#'
#' # base36 encoding of a large number
#' itostr(1e7, 36)
itostr = function(x, base = 10L) {
  x = asInteger(x, any.missing = FALSE, lower = 0L)
  base = asInt(base, na.ok = FALSE, lower = 2L, upper = 36L)
  .Call("itostr", x, base, PACKAGE = "BBmisc")
}
