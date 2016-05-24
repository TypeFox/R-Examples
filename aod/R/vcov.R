## Define a new generic vcov function if necessary
if(!isGeneric("vcov"))
  setGeneric(name = "vcov", def = function(object, ...) standardGeneric("vcov"))

setOldClass("geeglm")
setOldClass("geese")

## vcov method for objects of class "glimML"
setMethod(f = "vcov", signature = "glimML", definition = function(object, ...) {
  coeff <- coef(object)
  nc <- length(coeff)
  vc <- as.matrix(object@varparam[seq(nc), seq(nc)])
  cn <- names(coeff)
  dimnames(vc) <- list(cn, cn)
  vc
  })

## vcov method for objects of class "glimQL"
setMethod(f = "vcov", signature = "glimQL", definition = function(object, ...) vcov(object@fm))

# vcov method for objects of class geeglm
setMethod(f = "vcov", signature = "geeglm", definition = function(object, ...){
  b <- coef(object)
  vb <- object$geese$vbeta
  if(length(b) != NCOL(vb))
    stop("The number of coefficients does not match the dimension of the var-cov matrix.")
  nam <- names(b)
  dimnames(vb) <- list(nam, nam)
  vb
  })

#vcov.geeglm <- function(object, ...){
#  b <- coef(object)
#  vb <- object$geese$vbeta
#  if(length(b) != NCOL(vb))
#    stop("The number of coefficients does not match the dimension of the var-cov matrix.")
#  nam <- names(b)
#  dimnames(vb) <- list(nam, nam)
#  vb
#  } 

## vcov method for objects of class geese
setMethod(f = "vcov", signature = "geese", definition = function(object, ...){
  b <- object$beta
  vb <- object$vbeta
  if(length(b) != NCOL(vb))
    stop("The number of coefficients does not match the dimension of the var-cov matrix.")
  nam <- names(b)
  dimnames(vb) <- list(nam, nam)
  vb
  })

#vcov.geese <- function(object, ...){
#  b <- object$beta
#  vb <- object$vbeta
#  if(length(b) != NCOL(vb))
#    stop("The number of coefficients does not match the dimension of the var-cov matrix.")
#  nam <- names(b)
#  dimnames(vb) <- list(nam, nam)
#  vb
#  } 
