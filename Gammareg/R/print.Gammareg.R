print.Gammareg <-
function(x, ...)
{
  
  cat (" \n            ################################################################
            ###                  Classic Gamma Regression                ###
            ################################################################ \n")
  
  cat("\n Call: \n")
  
  print(x$call)
  cat("\n Coefficients: \n")
  print(x$coefficients)
  
}
