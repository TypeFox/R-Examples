print.lewbel.model <-
function(x, ...){
  cat("Call:\n")
  print(x$call)
  cat("\nCoefficients and Standard Errors:\n")
  printCoefmat(x$coef.est, P.values=T, has.Pvalue=T)
  cat("\nOveridentification Test:\n")
  print(x$j.test)
  cat("\nPartial F-test Statistics for Weak IV Detection:\n")
  print(x$f.test.stats)
  cat("\nNumber of Observations:", x$num.obs, "\n")
}
