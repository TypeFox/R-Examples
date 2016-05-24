print.Mort2Dsmooth <-
function(x, digits=max(3, getOption("digits")-3), ...){
  if(!is.null(cl <- x$call)){
    cat("Call:\n")
    dput(cl, control=NULL)
  }
  cat("\nNumber of Observations              :", x$n*x$m, "\n")
  cat("                    of which x-axis :", x$m, "\n")
  cat("                             y-axis :", x$n, "\n")
  cat("Effective dimension                 :", format(round(x$df, 4)), "\n")
  cat("(Selected) smoothing parameters      ", "\n")
  cat("                         over x-axis:", format(signif(x$lambdas[1], 5)), "\n")
  cat("                         over y-axis:", format(signif(x$lambdas[2], 5)), "\n")
  cat("Bayesian Information Criterion (BIC):", signif(x$bic, 5))
  cat("\n")
  invisible(x)
}
