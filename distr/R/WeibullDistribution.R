
################################
##
## Class: WeibullParameter
##
################################

## Access Methods
setMethod("shape", "WeibullParameter", 
           function(object) object@shape)
setMethod("scale", "WeibullParameter", 
           function(x, center = TRUE, scale = TRUE) x@scale)
           ### odd arg-list due to existing function in base package 

## Replace Methods
setReplaceMethod("shape", "WeibullParameter", 
                  function(object, value){ object@shape <- value; object})
setReplaceMethod("scale", "WeibullParameter", 
                  function(object, value){ object@scale <- value; object})

setValidity("WeibullParameter", function(object){
  if(length(shape(object)) != 1)
    stop("shape has to be a numeric of length 1")    
  if(shape(object) <= 0)
    stop("shape has to be positive")
  if(length(scale(object)) != 1)
    stop("scale has to be a numeric of length 1")      
  if(scale(object) <= 0)
    stop("scale has to be positive")
  else return(TRUE)
})

################################
##
## Class: Weibull distribution
##
################################

Weibull <- function(shape = 1, scale = 1) 
                    new("Weibull", shape = shape, scale = scale)

## wrapped access methods
setMethod("shape", "Weibull", 
           function(object) shape(param(object)))
setMethod("scale", "Weibull", 
           function(x, center = TRUE, scale = TRUE) scale(param(x)))
           ### odd arg-list due to existing function in base package 

## wrapped replace methods
setMethod("shape<-", "Weibull", 
           function(object, value) 
                    new("Weibull", shape = value, scale = scale(object)))
setMethod("scale<-", "Weibull", 
           function(object, value) 
                    new("Weibull", shape = shape(object), scale = value))

setMethod("*", c("Weibull","numeric"),
          function(e1, e2){
                if(isTRUE(all.equal(e2,0)))
                   return(new("Dirac", location = 0, .withArith = TRUE))
                if(e2 > 0)
                   return(new("Weibull", shape = shape(e1),
                               scale = scale(e1) * e2, .withArith = TRUE))
                return(-1 * as(Weibull(shape = shape(e1),
                                    scale = scale(e1) * (-e2)),
                               "AbscontDistribution")
                      )
          })
