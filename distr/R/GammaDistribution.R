
################################
##
## Class: GammaParameter
##
################################


## Access Methods
setMethod("shape", "GammaParameter", function(object) object@shape)
setMethod("scale", "GammaParameter", 
           function(x, center = TRUE, scale = TRUE) x@scale)
           ### odd arg-list due to existing function in base package 

## Replace Methods
setReplaceMethod("shape", "GammaParameter", 
                  function(object, value){ object@shape <- value; object})
setReplaceMethod("scale", "GammaParameter", 
                  function(object, value){ object@scale <- value; object})


validGammaParameter <- function(object){
  if(length(shape(object)) != 1)
    stop("shape has to be a numeric of length 1")    
  if(shape(object) <= 0)
    stop("shape has to be positive")
  if(length(scale(object)) != 1)
    stop("scale has to be a numeric of length 1")    
  if(scale(object) <= 0)
    stop("scale has to be positive")
  else return(TRUE)
}

setValidity("GammaParameter", validGammaParameter)


################################
##
## Class: gamma distribution
##
################################

Gammad <- function(shape = 1, scale = 1) 
          new("Gammad", shape = shape, scale = scale)

## wrapped access methods
setMethod("shape", "Gammad", function(object) shape(param(object)))
setMethod("scale", "Gammad", 
           function(x, center = TRUE, scale = TRUE) scale(param(x)))
           ### odd arg-list due to existing function in base package 

## wrapped replace methods
setMethod("shape<-", "Gammad", 
           function(object, value) 
                    new("Gammad", shape = value, scale = scale(object)))
setMethod("scale<-", "Gammad", 
           function(object, value) 
                    new("Gammad", shape = shape(object), scale = value))

