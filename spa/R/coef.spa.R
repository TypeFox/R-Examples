"coef.spa" <-function(object,...){
  if(class(object)!="spa"){
    stop("Error:  x is not of type spa")
  }
  return(object$model$xstr$coef[,1])
}

