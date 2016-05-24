print.glinternet.cv = function(x, ...){

  #print the call
  cat("Call: ")
  print(x$call)

  #print best lambda, minimum cv error
  cat("Results of ", x$nFolds, "-fold cross validation:\n")
  cat("Minimum CV error of ", min(x$cvErr), " at lambda = ", x$lambda[which.min(x$cvErr)[1]], "\n")
  
  #print the cv errors and standard errors
  output = data.frame(lambda=signif(x$lambda, 3), cvErr=x$cvErr, cvErrStd=x$cvErrStd)
  print(output)
}
