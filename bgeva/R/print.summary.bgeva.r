print.summary.bgeva <- function(x,digits = max(3, getOption("digits") - 3),
                                             signif.stars = getOption("show.signif.stars"),...){

  cat("\nFamily: BGEVA \nEquation: ") 
  print(x$formula)
  cat("\n") 
  cat("Parametric coefficients:\n") 
  printCoefmat(x$tableP,digits = digits, signif.stars = signif.stars,na.print = "NA",...)
  cat("\n")

    if(x$l.sp!=0){
    cat("Smooth components' approximate significance:\n")
    printCoefmat(x$tableNP,digits = digits, signif.stars = signif.stars,has.Pvalue = TRUE,na.print = "NA",cs.ind = 1,...)
    cat("\n")
    }

    cat("n = ",x$n,"  tau = ",format(x$tau,digits=3),"  total edf = ",format(x$t.edf,digits=3),"\n\n", sep="") 
              
invisible(x)
                
}
