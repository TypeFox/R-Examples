#' Prints object to a string / character vector.
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
#' x = data.frame(a = 1:2, b = 3:4)
#' str(printToChar(x))
printToChar = function(x, collapse = "\n") {
  rval = NULL
  con = textConnection("rval", "w", local = TRUE)
  sink(con)
  on.exit({
    sink()
    close(con)
  })
  print(x)
  if (!is.null(collapse))
    paste(rval, collapse = collapse)
  else
    rval
}
