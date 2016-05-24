#' hrr provides a simple function, \code{\link{hr}}, to compose and print beautiful horizontal rules.
#'
#' @section Package options:
#'
#' It uses the following \code{options} to configure behaviour:
#' \itemize{
#'   \item \code{hrr.symbol}: default character symbol used by \code{\link{hr}} function
#' }
#' @docType package
#' @name hrr
NULL

.onLoad <- function(libname, pkgname) {
  ops <- options()
  hrr_ops <- list(
    hrr.symbol = '#'
  )
  to_set <- !(names(hrr_ops) %in% names(ops))
  if (any(to_set)) options(hrr_ops[to_set])
  invisible()
}
