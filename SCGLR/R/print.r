#' @export
#' @title Print SCGLR object
#' @description Prints inertia per component and deviance for each Y.
#' @method print SCGLR
#' @param x object of class 'SCGLR', usually a result of running \code{\link{scglr}}.
#' @param \dots Not used.
print.SCGLR <- function(x, ...) {
  cat("\nCall:  ", paste(deparse(x$call), sep = "\n", collapse = "\n"), sep = "","\n")
  cat("\nInertia:\n")
  print.default(x$inertia,print.gap=2)
  cat("\nDeviance:\n")
  print.default(x$deviance,print.gap=2)
  cat("\n")
  invisible(x)
}
