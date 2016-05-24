print.summary.Mort2Dsmooth <-
function(x, digits=max(3, getOption("digits")-3), ...){
  if(!is.null(cl <- x$call)){
    cat("Call:\n")
    dput(cl, control=NULL)
  }
  cat("\nNumber of Observations                      :", x$n*x$m, "\n")
  cat("                            of which x-axis :", x$m, "\n")
  cat("                                     y-axis :", x$n, "\n")
  cat("Effective dimension                         :", format(signif(x$df, 4)), "\n")
  cat("(Selected) smoothing parameters              ", "\n")
  cat("                                 over x-axis:", format(signif(x$lambdas[1], 5)), "\n")
  cat("                                 over y-axis:", format(signif(x$lambdas[2], 5)), "\n")
  cat("Bayesian Information Criterion (BIC)        :", signif(x$bic,5), "\n")
  cat("Akaike's Information Criterion (AIC)        :", signif(x$aic,5), "\n")
  cat("(Estimated) overdispersion parameter (psi^2):", signif(x$psi2,3), "\n")
  cat("\nResiduals:\n", sep="")
  nam <- c("Min", "1Q", "Median", "3Q", "Max")
  rq <- structure(quantile(x$residuals), names = nam)
  print(rq, digits = digits, ...)
  cat("\nSettings and control:\n")
  cat("  x-axis:", x$ndx[1] + x$deg[1], "B-splines of degree", x$deg[1], "- differences of order", x$pord[1], "\n")
  cat("  y-axis:", x$ndx[2] + x$deg[2], "B-splines of degree", x$deg[2], "- differences of order", x$pord[2], "\n")
  cat("  convergence tolerance :", format(signif(x$tolerance, 5)))
  cat("\n")
  invisible(x)
}
