"dim.spa" <-function(x){
  if(class(x)!="spa"){
    stop("Error:  x is not of type spa")
  }
  dims=x$model$dims
  return(dims[-3])
}

