print.summary.Bayesianbetareg <-
function(x, ...){
  
  cat (" \n            ################################################################
            ###                  Bayesian Beta Regression                ###
            ################################################################ \n")
  
  cat("\n Call: \n")
  print(x$call)
  cat("\n")
  
  printCoefmat(x$coefficients, digits=4)
  
  cat("\n Deviance: \n")
  print(x$Deviance)
  
  cat("\n AIC: \n")
  print(x$AIC)
  
  cat("\n BIC: \n")
  print(x$BIC)
}
