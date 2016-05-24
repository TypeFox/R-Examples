print.stratamean <-
function(x,...)
{
   cat("\nstratamean object: Stratified sample mean estimate\n", sep="")
   if(x$call$fpc == TRUE) cat("With ")
   else cat("Without ")
   cat("finite population correction.")
   cat("\nMean estimate: ", round(x$mean,4), "\n", sep="")
   cat("Standard error: ", round(x$se,4), "\n", sep="")
   cat(100 * x$call$level,"% confidence interval: [",round(x$ci[1],4),",",round(x$ci[2],4),"]\n\n", sep="")
   invisible(x)
}
