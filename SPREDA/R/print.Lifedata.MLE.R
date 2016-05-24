print.Lifedata.MLE <-
function(x,...){
  m=length(x$coef)
  coef=c(x$coef[1:(m-1)], exp(x$coef[m]))
  names(x)[m]="sigma"
  cat("Call:\n")
  print(x$call)
  cat("\n")
  cat("Coefficients:\n")
  print(coef)
  cat("\n")
  cat("Loglikelihod:\n")
  cat(as.numeric(-x$min), " (df=", sum(coef!=0), ")", sep="")
}
