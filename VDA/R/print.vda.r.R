print.vda.r <-
function (x, ...)
{
  cat ("\n Call: \n")
  print (x$call)

  cat ("# of cases=", x$cases, "\n")
  cat ("# of classes=", x$classes, "\n")
  cat ("# of features=", x$features, "\n")
  
  cat ("\n Lambda used:", x$lambda, "\n")
  
  cat ("\n Estimated Coefficients: \n")
  print (x$coefficient)  
  
  cat ("\n Predicted classification: \n")
  print (x$predicted)
  
  cat ("\n Training Error: \n")
  print (x$training_error_rate)
}
