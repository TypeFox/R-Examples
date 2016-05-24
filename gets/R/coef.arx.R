coef.arx <-
function(object, spec=NULL, ...)
{
  ##spec argument:
  if(is.null(spec)){
    spec <- switch(as.character(object$call)[1],
      arx="both", getsm="mean", getsv="variance")
    #spec <- "both"
  }else{
    spec.type <- c("mean", "variance", "both")
    which.type <- charmatch(spec, spec.type)
    spec <- spec.type[which.type]
  } #end if(..)else(..)

  ##mean:
  result1 <- NULL
  if(spec=="mean" || spec=="both"){
    if(!is.null(object$mean.results)){
      result1 <- object$mean.results[,1]
      names(result1) <- rownames(object$mean.results)
    }
  } #end if(spec==..)

  ##variance:
  result2 <- NULL
  if(spec=="variance" || spec=="both"){
    if(!is.null(object$variance.results)){
      result2 <- object$variance.results[,1]
      names(result2) <- rownames(object$variance.results)
      if(!is.null(object$Elnz2)){
        result2 <- c(result2,object$Elnz2)
        names(result2)[length(result2)] <- "Elnz2"
      }
    } #end if(..)else(..)
  } #end if(spec==..)

  result <- c(result1,result2)
  return(result)
}
