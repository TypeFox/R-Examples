print.logitor <-
function(x, ...){
  cat("Call:\n")
  print(x$call)
  cat("\nOdds Ratio:\n")
  printCoefmat(x$oddsratio,P.values=T,has.Pvalue=T)
}
