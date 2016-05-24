summary.l1.reg <-
function (object, ...)
{
  cat ("\n Call: \n")
  print (object$call)
  cat ("\n")
  
  cat ("\n Feature Matrix \n")
  print (object$X)
  cat ("\n")

  cat ("\n Outcome for Cases \n")
  print (object$Y)
  cat ("\n")
  
  cat ("\n Residuals \n")
  print (object$residual)
  cat ("\n")
  
  cat ("\n Sum of Residuals \n")
  print (object$L1)
  cat ("\n")
  
  cat ("# of cases=", object$cases, "\n")
  cat ("# of predictors=", object$predictors, "\n")
  
  cat ("\n Lambda used:", object$lambda, "\n")
  cat ("\n")
  
  cat ("\n Intercept: \n")
  print (object$intercept)
  
  cat ("\n Esimtated Coefficients: \n")
  print (object$estimate)
  cat ("\n")
  
  cat ("\n Number of Active Variables: \n")
  print (object$nonzeros)
  cat ("\n")
    
  cat ("\n Selected Variables with Nonzero Coefficients: \n")
  print (rownames(object$X)[object$selected+1])
  }
