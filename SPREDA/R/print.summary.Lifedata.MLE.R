print.summary.Lifedata.MLE <-
function(x,...){
  cat("Call:\n")
  print(x$call)
  cat("\n")
  cat("Parameters:\n")
  print(x$coefmat)
  cat("\n")
  cat("Loglikelihod:\n")
  print(as.numeric(-x$min))
  cat("\n")
  cat("Covariance matrix:\n")
  print(x$vcov)
  cat("\n")
  cat("Survival probability:\n")
  print(x$surv)
}
