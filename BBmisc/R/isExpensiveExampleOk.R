#' Conditional checking for expensive examples.
#'
#' Queries environment variable \dQuote{R_EXPENSIVE_EXAMPLE_OK}.
#' Returns \code{TRUE} iff set exactly to \dQuote{TRUE}.
#' This allows conditional checking of expensive examples in packages
#' via R CMD CHECK, so they are not run on CRAN, but at least
#' on your local computer.
#' A better option than \dQuote{dont_run} in many cases, where such examples
#' are not checked at all.
#'
#' @return [\code{logical(1)}].
#' @export
#' @examples
#' # extremely costly random number generation, that we dont want checked on CRAN
#' if (isExpensiveExampleOk()) {
#'   runif(1)
#' }
isExpensiveExampleOk = function() {
  Sys.getenv("R_EXPENSIVE_EXAMPLE_OK") == "TRUE"
}
