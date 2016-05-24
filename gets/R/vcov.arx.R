vcov.arx <-
function(object, spec=NULL, ...)
{

  ##spec argument:
  if(is.null(spec)){
    if(!is.null(object$mean.results)){
      spec <- "mean"
    }
    if(is.null(object$mean.results)
      && !is.null(object$variance.results) ){
      spec <- "variance"
    }
  }else{
    spec.type <- c("mean", "variance")
    which.type <- charmatch(spec, spec.type)
    spec <- spec.type[which.type]
  }

  ##mean:
  if(spec=="mean"){
    result <- object$vcov.mean
  }

  ##variance:
  if(spec=="variance"){
    result <- object$vcov.var
  }

#  ##check and change if 0 x 0?:
#  if(all(dim(result)==0)){ result <- NULL }

  return(result)
}
