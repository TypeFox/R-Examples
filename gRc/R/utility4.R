
getSlot<-function(object, slot){
  object[[slot]]
}

getSlot<-function(object, slot){
  if(is.null(slot))
    return(object)
  
  return(object[[slot]])
}

dataRep <- function(object,slot=NULL){
  if (is.null(slot))
    return(getSlot(object,"dataRep"))
  getSlot(object,"dataRep")[[slot]]
}

intRep <- function(object,slot=NULL){
  if (is.null(slot))
    return(getSlot(object,"intRep"))
  getSlot(object,"intRep")[[slot]]
}

fitInfo <- function(object,slot=NULL){
  if (is.null(slot))
    return(getSlot(object,"fitInfo"))
  getSlot(object,"fitInfo")[[slot]]
}

coef.rcox <- function(object, ...){
  co  <- fitInfo(object,"coef")
  co
}

getcc <- function(object,type){
  if (missing(type))
    list(vcc=object$vcc, ecc=object$ecc)
  else {
    switch(type,
           "ecc"={object$ecc},
           "vcc"={object$vcc})
  }
}

getecc <- function(object){
  object$ecc
}

getvcc <- function(object){
  object$vcc
}

getedges <- function(object,complement=FALSE){
  ans <- ecc2edges(getecc(object))
  if (complement){
    eAll <- names2pairs(getSlot(object,"nodes"))
    ans  <- setdiffLL(eAll, ans)
  }
  ans
}


print.colourClass <- function(x,...){
  xf <- names2formula(x)
  xs <- formula2string(xf)
  mapply(function(n,xxx) cat(n,xxx,"\n"), names(xs),xs)
  return(invisible(x))
}





tocc <- function(v){
  if(length(v)==0)
    return(NULL)
  as.cclist(
  lapply(v, function(x) {
    if (length(x)==1)
      as.cc(as.atom(x))
    else
      as.cc(lapply(x, as.atom))
    })
  )
}


cc2str <- function(cc){
  paste(sapply(cc, toLisp),collapse='')  
}

.addccnames <- function(x, type){
  if (length(x)){
    names(x) <- paste(type,paste(1:length(x)),sep="")
    class(x) <- c("colourClass","list")
    x
  } else {
    NULL
  }
}

## ellK <- function(K, S, n){
##   value <- (n/2)*(log(det(K)) - sum(rowSums(K*S)))

##   ##diag(crossprod(K,S))))
##   return(value)
## }



dimension  <- function(m){
  length(c(getSlot(m,'vcc'), getSlot(m,'ecc')))
}

logL <- function(m){
  getSlot(m,'fitInfo')$logL
}

cholSolve <- function(ma)
  chol2inv(chol(  ma  ))
