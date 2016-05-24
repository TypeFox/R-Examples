print.summary.Gammareg <-
function(x, ...){
  
  cat (" \n            ################################################################
            ###                  Classic Gamma Regression                ###
            ################################################################ \n")
  
  cat("\n Call: \n")
  print(x$call)
  cat("\n")
  
  printCoefmat(x$coefficients, digits=4)
  
  cat("\n Covariance Matrix for Beta: \n")
  print(x$covB)
  
  cat("\n Covariance Matrix for Gamma: \n")
  print(x$covG)

  cat("\n AIC: \n")
  print(x$AIC)
  
  
  cat("\n Iteration: \n")
  print(x$iteration)
  

  cat("\n Convergence: \n")
  print(x$convergence)
}
