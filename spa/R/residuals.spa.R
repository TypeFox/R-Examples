"residuals.spa" <-function(object,verbose=FALSE,...){
  if(class(object)!="spa"){
    stop("Error:  x is not of type spa")
  }
  if(verbose)
    cat("The are for  are given for y[L]")
  return(object$model$res)
}

