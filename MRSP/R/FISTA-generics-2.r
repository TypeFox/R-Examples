################################################################################
#####    Generic functions for FISTA                                       #####
################################################################################
#####    Author: Wolfgang Pößnecker                                        #####
#####    Last modified: 17.01.2014, 16:31                                  #####
################################################################################

setGeneric("fistaProximal",
           function(coef, tuning, penweights, penindex, grpindex, Proximal.args, Proximal.control, ...)
           standardGeneric("fistaProximal"))
           
setGeneric("penalty",
           function(coef, tuning, penweights, penindex, grpindex, Proximal.args, ...)
           standardGeneric("penalty"))

setGeneric("lossapprox",
           function(lS, grad, coefdiff, Proximal.args, ...)
           standardGeneric("lossapprox"))

## standard method for lossapprox:           
setMethod("lossapprox", signature(coefdiff = "list"),
          function(lS, grad, coefdiff, Proximal.args, ...)
{
 lS + Reduce('+', Map(function(u,v){sum(u*v)}, grad, coefdiff))
}) 

setGeneric("proximterm",
           function(coefdiff, Proximal.args, ...)
           standardGeneric("proximterm"))

## standard method for proximterm:           
setMethod("proximterm", signature(coefdiff = "list"),
          function(coefdiff, Proximal.args, ...)
{
 Reduce('+', lapply(coefdiff, function(u) sum(u^2)))
})

setGeneric("useSPG",
           function(coef, penweights, Proximal.args, ...)
           standardGeneric("useSPG"))
           
setMethod("useSPG", signature(coef = "list"),
          function(coef, penweights, Proximal.args, ...){F})

setGeneric("SPGsmoothpen",
           function(coef, dualopt, penweights, grpindex, tuning, Proximal.args, ...)
           standardGeneric("SPGsmoothpen"))

setMethod("SPGsmoothpen", signature(coef = "list"),
          function(coef, dualopt, penweights, grpindex, tuning, Proximal.args, ...){0})

setGeneric("SPG",
           function(coef, tuning, penweights, grpindex, Proximal.args, doGrad, ...)
           standardGeneric("SPG"))
           
## standard method for SPG:
setMethod("SPG", signature(coef = "list"),
          function(coef, tuning, penweights, grpindex, Proximal.args, doGrad=T, ...)
{
 SPGlist <- list()
 if(doGrad){
  SPGlist$grad <- lapply(coef, function(u){if(is.matrix(u)){u <- matrix(0, nrow=nrow(u), ncol=ncol(u))}
                                           else if(is.vector(u)){u <- rep(0, length(u))}
                                           else stop("object coef had an unusable form in function 'SPG'")})
 }
 SPGlist$dualopt <- list()
 return(SPGlist)
})
           
setGeneric("updateEta",
           function(dat, coef, offset, weights, ...)
           standardGeneric("updateEta"))

setGeneric("updatePenweights",
           function(penweights, coef, weights, penindex, grpindex, Proximal.args, ...)
           standardGeneric("updatePenweights"))
           
setGeneric("logLik")

setGeneric("extract",
           function(object, what, ...)
           standardGeneric("extract"))  
