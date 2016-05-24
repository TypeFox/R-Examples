print.summary.factorQR <- function(x, digits = max(3, getOption("digits") - 3), ...) {

  cat("\nCall:\n")
  cat(paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n\n",
      sep = "") 

  cat("\n Lambda components:\n")
  printCoefmat(x$coef.table, digits = digits, na.print = "NA", ...)
  
  cat("\n Other components:\n")
  printCoefmat(x$cov.table, digits = digits, na.print = "NA", ...)
  
  cat("\nNumber of observations:", x$n.obs)
  cat("\nNumber of estimated parameters:", x$n.param)
  cat("\nNumber of stored MCMC draws:", x$n.draws)
  cat("\n\n")
  invisible(x)
}
