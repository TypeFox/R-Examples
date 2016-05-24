##
## Methods for returning the convergence codes of the component GARCH models
## =========================================================================
##
## Method definition for objects of class "GoGARCH"
##
setMethod(f = "converged", signature(object = "GoGARCH"), definition = function(object, ...){
  conv <- c(unlist(lapply(object@models, function(x) x@fit$convergence)))
  cnames <- paste("y", seq(along.with = conv), sep = "")
  names(conv) <- cnames
  return(conv)
})
##
## Method definition for objects of class "Goestica"
## "Goestica" extends directly "GoGARCH"
##
setMethod(f = "converged", signature(object = "Goestica"), definition = function(object){
  converged(as(object, "GoGARCH"))
})
##
## Method definition for objects of class "Goestmm"
## "Goestmm" extends directly "GoGARCH"
##
setMethod(f = "converged", signature(object = "Goestmm"), definition = function(object){
  converged(as(object, "GoGARCH"))
})
##
## Method definition for objects of class "Goestnls"
## "Goestnls" extends directly "GoGARCH"
##
setMethod(f = "converged", signature(object = "Goestnls"), definition = function(object){
  converged(as(object, "GoGARCH"))
})
##
## Method definition for objects of class "Goestml"
## "Goestml" extends directly "GoGARCH"
##
setMethod(f = "converged", signature(object = "Goestml"), definition = function(object){
  converged(as(object, "GoGARCH"))
})

