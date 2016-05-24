#' Print Method for summary.mpm Objects
#' @param x object of class summary.mpm
#' @param digits minimum number of significant digits to print, defaults to 2
#' @param what one of \code{"columns"} (default), \code{"rows"} or \code{"all"}, specifying respectively
#'   whether columns, rows or both need to be printed
#' @param ... further arguments for the print method
#' @return x is returned invisibly
#' @seealso  \code{\link{print.default}}
#' @S3method print summary.mpm
#' @method print summary.mpm
#' @export
print.summary.mpm <- function(x, digits = 2, 
    what = c("columns", "rows", "all"), ...){
  if (missing(x)) stop("Argument \"x\" is missing, with no default")
  what <- match.arg(what)
  cat("\nCall:\n")
  dput(x$call)
  cat("\n", length(x$Columns$Posit), " columns and ", length(x$Row$Posit), " rows.\n")
  cat("\n", sum(x$Columns$Posit), " columns and ", sum(x$Rows$Posit), " rows positioned.\n")
  cat("\nContributions:\n")
  cx <- format(round(x$VPF,digits), digits = digits)
  print(cx, quote = FALSE, ...)
  if (what == "columns" || what == "all"){
    cat("\nColumns:\n")
    cx <- format(round(x$Columns, digits), digits = digits)
    print(cx, quote = FALSE, ...)
  }
  if (what == "rows" || what == "all"){
    cat("\nRows:\n")
    cx <- format(round(x$Rows, digits), digits = digits)
    print(cx, quote = FALSE, ...)
  }
  invisible(x)
}
