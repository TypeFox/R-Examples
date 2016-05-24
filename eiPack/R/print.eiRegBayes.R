print.eiRegBayes <- function(x, digits = max(2, getOption("digits") - 4), ...) {
  cat("\nFormula: ", deparse(x$call$formula), "\n")
  cat("Total sims: ", dim(x$draws)[3], "\n\n")
  cat("Estimated internal cells: (averaged across simulations)\n")
  print.default(format(apply(x$draws, c(1,2), mean), digits = digits), 
                print.gap = 2, quote = FALSE, ...)
  invisible(x)
}
