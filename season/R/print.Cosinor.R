## print.Cosinor.R
## Prints basic results from Cosinor

print.Cosinor<-function(x, ...){

  ## Checks
  if (class(x)!="Cosinor"){stop("Object must be of class 'Cosinor'")} 

  ## Use GLM function ###
  print(x$glm, ...)
} # end of function
