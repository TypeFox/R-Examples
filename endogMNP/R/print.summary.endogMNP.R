print.summary.endogMNP <- function(x, digits = max(3, getOption("digits") - 3), ...) {

  cat("\nCall:\n")
  cat(paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n\n",
      sep = "") 

  cat("\nCoefficients:\n")
  printCoefmat(x$coef.table, digits = digits, na.print = "NA", ...)
  
  cat("\nCovariances:\n")
  printCoefmat(x$cov.table, digits = digits, na.print = "NA", ...)
  
  cat("\nSelection base category:", x$selBase)  
	cat("\nOutcome base category:", x$outBase)  
  cat("\nNumber of observations:", x$n.obs)
  cat("\nNumber of estimated parameters:", x$n.param)
  cat("\nNumber of stored MCMC draws:", x$n.draws)
  cat("\n\n")
  invisible(x)
}
