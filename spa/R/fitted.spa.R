"fitted.spa" <-function(object,verbose=FALSE,...){
  if(class(object)!="spa"){
    stop("Error:  x is not of type spa")
  }
  if(verbose)
    cat("The fitted values are given for", dim(object)[1]," values that which", dim(object)[2]," are labeled\nIn addition, they are ordered as given by the inputed response\n")
  return(object$model$fit)
}

