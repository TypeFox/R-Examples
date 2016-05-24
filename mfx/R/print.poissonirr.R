print.poissonirr <-
function(x, ...){
  cat("Call:\n")
  print(x$call)
  cat("\nIncidence-Rate Ratio:\n")
  printCoefmat(x$irr, P.values=T,has.Pvalue=T)
}
