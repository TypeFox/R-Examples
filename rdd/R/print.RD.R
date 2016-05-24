#'Print the Regression Discontinuity
#' 
#' Print a very basic summary of the regression discontinuity
#' 
#' @method print RD
#' @param x \code{rd} object, typically the result of \code{\link{RDestimate}}
#' @param digits number of digits to print
#' @param ... unused
#' @include RDestimate.R
#' @export
#' @author Drew Dimmery <\email{drewd@@nyu.edu}>
print.RD<-function(x,digits=max(3, getOption("digits") - 3),...) {
  cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"), 
      "\n\n", sep = "")
  
  cat("Coefficients:\n")
  print.default(format(x$est,digits = digits),print.gap=2,quote=FALSE)
  cat("\n")
  invisible(x)
  
}