print.eiRegBayesSum <- function(x, digits = max(2, getOption("digits") - 4), ...) {
  cat("\nFormula: ", deparse(x$call$formula), "\n")
  cat("Total sims: ", x$sims, "\n\n")
  cat("Estimated internal cells: (across simulations)\n")
  print.default(format(x$coef, digits = digits), 
                print.gap = 2, quote = FALSE, ...)
  invisible(x)
}
