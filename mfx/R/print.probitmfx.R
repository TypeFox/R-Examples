print.probitmfx <-
function(x, ...){
  cat("Call:\n")
  print(x$call)
  cat("\nMarginal Effects:\n")
  printCoefmat(x$mfxest,P.values=T,has.Pvalue=T)
  if(length(x$dcvar)!=0){
    cat("\ndF/dx is for discrete change for the following variables:\n\n")
    print(x$dcvar)
  }
}
