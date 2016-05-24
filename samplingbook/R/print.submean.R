print.submean <-
function(x,...)
{
  cat("\nsubmean object: Sub-sample mean estimate\n", sep="")
  if(x$call$fpc == FALSE) cat("Without ")
  else cat("With ")
  cat("finite population correction.") 
  cat("\nUsing method: ", x$call$method, "\n", sep="")
  cat("\nMean estimate: ", round(x$mean,4), "\n", sep="")
  cat("Standard error: ", round(x$se,4), "\n", sep="")
  cat(100 * x$call$level,"% confidence interval: [",round(x$ci[1],4),",",round(x$ci[2],4),"]\n\n", sep="")
  invisible(x)
}
