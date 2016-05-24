print.deglmx <-
function(x, ...){
  
  cat("Degradation fit:\n")
  print(x$type)
  cat("\n")
  
  cat("Degradation trend:\n")
  print(x$ytrend)
  cat("\n")
  
  cat("Dynamic covariates:\n")
  print(x$dyncovnames)
  cat("\n")
  
  cat("Convergent:\n")
  print(x$fit$conv)
  cat("\n")
  
  cat("Coefficients:\n")
  print(x$fit$coef)
  cat("\n")
  
  cat("Loglikelihood:\n")
  print(x$fit$loglik) 
  cat("\n")
  
  cat("Random effects:\n")
  print(x$fit$ran.eff)
}
