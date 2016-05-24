fitted.arx <-
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
    spec.type <- c("mean", "variance", "both")
    which.type <- charmatch(spec, spec.type)
    spec <- spec.type[which.type]
  }

  result <- NULL

  ##mean:
  if(spec=="mean"){
    result <- object$mean.fit
  }

  ##variance:
  if(spec=="variance"){
    result <- object$var.fit
  }

  ##both:
  if(spec=="both"){
    if(!is.null(object$mean.results)
      && !is.null(object$variance.results) ){
      result <- cbind(object$mean.fit, object$var.fit)
      colnames(result) <- c("yhat","sigma2hat")
    }
  }

  return(result)
}
