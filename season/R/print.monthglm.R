## print.monthglm.R
## Prints basic results from monthglm
print.monthglm<-function(x, ...){
  ## Checks
  if (class(x)!="monthglm"){stop("Object must be of class 'monthglm'")} 
  ## Use GLM function ###
  print(x$glm, ...)
} # end of function

