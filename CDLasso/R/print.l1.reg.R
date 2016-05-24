print.l1.reg <-
function (x, ...)
{
  cat ("\n Call: \n")
  print (x$call)
  
  cat ("# of cases=", x$cases, "\n")
  cat ("# of predictors=", x$predictors, "\n")
  
  cat ("\n Lambda used:", x$lambda, "\n")
  cat ("\n")
  
  cat ("\n Intercept: \n")
  print (x$intercept)
  
  cat ("\n Selected Coefficient Estimates: \n")
  print (cbind(Predictor=rownames(x$X)[x$selected+1],Estimate=x$estimate[x$selected]),justify="centre")
  
  cat ("\n Number of Active Variables: \n")
  print (x$nonzeros)
}
