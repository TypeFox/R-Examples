print.Bayesianbetareg <-
function(x, ...)
{
  
  cat (" \n            ################################################################
            ###                  Bayesian Beta Regression                ###
            ################################################################ \n")
  
  cat("\n Call: \n")
  
  print(x$call)
  cat("\n Coefficients: \n")
  print(x$coefficients)
  
}
