print.summary.Mort1Dsmooth <-
function(x, digits=max(3, getOption("digits")-3), ...){
  if(!is.null(cl <- x$call)){
    cat("Call:\n")
    dput(cl, control=NULL)
  }
  cat("\nNumber of Observations                  :", x$n, "\n")
  cat("Effective dimension                     :", format(signif(x$df, 4)), "\n")
  cat("(Selected) smoothing parameter          :", format(signif(x$lambda, 5)), "\n")
  cat("Bayesian Information Criterion (BIC)    :", signif(x$bic,5), "\n")
  cat("Akaike's Information Criterion (AIC)    :", signif(x$aic,5), "\n")
  cat("(Estimated) dispersion parameter (psi^2):", signif(x$psi2,3), "\n")
  cat("\nResiduals:\n", sep="")
  nam <- c("Min", "1Q", "Median", "3Q", "Max")
  rq <- structure(quantile(x$residuals), names = nam)
  print(rq, digits = digits, ...)
  cat("\nSettings and control:\n")
  cat("  number of B-splines    :", x$ndx + x$deg, "\n")
  cat("  degree of the B-splines:", x$deg, "\n")
  cat("  order of differences   :", x$pord, "\n")
  cat("  convergence tolerance  :", format(signif(x$tolerance, 5)))
  cat("\n")
  invisible(x)
}
