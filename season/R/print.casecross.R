## print.casecross.R
## Prints basic results from casecross
## Oct 2009
print.casecross<-function(x, ...){
## Check
  if (class(x) != "casecross"){
    stop("Object must be of class 'casecross'")
  } 
## Use print.coxph
  if (class(x$c.model) != "coxph"){
    stop("Conditional logistic regression model object 'c.model' must be of class 'coxph'")
  }    
  print(x$c.model, ...)
} # end of function
