summary.vda.le <-
function (object, ...)
{
  cat ("\n Call: \n")
  print (object$call)
  cat ("\n")
  
  cat ("# of cases=", object$cases, "\n")
  cat ("# of classes=", object$classes, "\n")
  cat ("# of features=", object$features, "\n")
  
  cat ("\n Lambdas used:", object$lambda, "\n")
  
  cat ("\n Predicted classification: \n")
  print (object$predicted)
  
  cat ("\n Training error: \n")
  print (object$training_error_rate)
  cat ("\n")
  
  cat ("\n Coefficients: \n")
  print (object$coefficient)
  cat ("\n")
  
  cat ("\n Number of Active Variables: \n")
  print (object$nonzeros)
  
  cat ("\n Selected Variables with Nonzero Coefficients: \n")
  print (names(object$feature)[object$selected+1])
  }
