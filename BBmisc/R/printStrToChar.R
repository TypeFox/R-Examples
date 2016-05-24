#' Print \code{str(x)} of an object to a string / character vector.
#'
#' @param x [any]\cr
#'   Object to print
#' @param collapse [\code{character(1)}]\cr
#'   Used to collapse multiple lines.
#'   \code{NULL} means no collapsing, vector is returned.
#'   Default is \dQuote{\\n}.
#' @return [\code{character}].
#' @export
#' @examples
#' printStrToChar(iris)
printStrToChar = function(x, collapse = "\n") {
  d = printToChar(str(x), collapse = NULL)
  # remove NULL from str
  collapse(d[-length(d)], collapse)
}
