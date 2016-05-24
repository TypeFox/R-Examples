# extract loadings from sempls object
plsLoadings <- function(object){
  UseMethod("plsLoadings", object)
}

plsLoadings.sempls <- function(object)
{
  plsLoadings <- object$cross_loadings
  # attribute for the measurement model
  attr(plsLoadings, "mm") <- object$model$M
  class(plsLoadings) <- c("plsLoadings")
  return(plsLoadings)
}

print.plsLoadings <- function(x, type=c("discriminant", "outer", "cross"),
                              cutoff=NULL, reldiff=0.2, na.print=".", digits=2,
                              abbreviate=FALSE,...)
{
  type <- match.arg(type)
  # to check discriminant validity
  if(type=="discriminant"){
    if(reldiff < 0 | reldiff > 1)
      stop("Argument 'reldiff' only accepts values in the intervall [0,1].")
    cross <- x
    outer <- cross
    outer[attr(cross, "mm")!=1] <- NA
    maxv <- apply(outer, 1, max, na.rm=TRUE)
    if(!is.null(cutoff)) cross[cross < cutoff] <- NA
    mind <- cross <= (maxv - reldiff * maxv)
    cross[mind] <- NA
    if(abbreviate) dimnames(cross) <- lapply(dimnames(cross), abbreviate, ...)
    print.table(cross, na.print=na.print, digits=digits, ...)
    invisible(x)
  }
  if(type=="outer"){
    cross <- x
    outer <- cross
    outer[attr(cross, "mm")!=1] <- NA
    if(!is.null(cutoff)){
      outer[outer <= cutoff] <- NA
      print.table(x=outer, na.print=na.print, digits=digits,...)
      invisible(x)
    }
    else{
      print.table(x=outer, na.print=na.print, digits=digits, ...)
      invisible(cross)
    }
  }
  if(type=="cross"){
    cross <- x
    if(!is.null(cutoff)) cross[cross <= cutoff] <- NA
    print.table(cross, na.print=na.print, digits=digits, ...)
    invisible(x)
  }
}

