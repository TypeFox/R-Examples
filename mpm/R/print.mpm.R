#' Print Method for mpm Objects
#' @param x object of class mpm
#' @param digits minimum number of significant digits to be printed
#' @param ... further arguments for the print method (for printing the contributions)
#' @seealso \code{\link{print.default}}
#' @return x is returned invisibly
#' @S3method print mpm
#' @method print mpm
#' @export
print.mpm <- function(x, digits = 3, ...){
  cat("Call:\n")
  dput(x$call)
  cat("\nContributions:\n")
  print(round(x$contrib, digits), ...)
  cat("\n", length(x$pos.column), " columns and ", length(x$pos.row), " rows.\n")
  cat("\n", sum(x$pos.column), " columns and ", sum(x$pos.row), " rows positioned.\n")
  invisible(x)
}
