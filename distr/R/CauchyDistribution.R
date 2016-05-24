
################################
##
## Class: CauchyParameter
##
################################


## Access methods
setMethod("location", "CauchyParameter", function(object) object@location)
setMethod("scale", "CauchyParameter", 
           function(x, center = TRUE, scale = TRUE) x@scale)
           ### odd arg-list due to existing function in base package 

## Replace Methods
setReplaceMethod("location", "CauchyParameter", 
                  function(object, value){ object@location <- value; object})
setReplaceMethod("scale", "CauchyParameter", 
                  function(object, value){ object@scale <- value; object})

setValidity("CauchyParameter", function(object){
  if(length(location(object)) != 1)
    stop("location has to be a numeric of length 1")    
  if(length(scale(object)) != 1)
    stop("scale has to be a numeric of length 1")    
  if(scale(object) <= 0)
    stop("scale has to be positive")
  else return(TRUE)
})



################################
##
## Class: CauchyDistribution
##
################################

Cauchy <- function(location = 0, scale = 1) 
                   new("Cauchy", location = location, scale = scale)

## wrapped access methods
setMethod("location", "Cauchy", function(object) location(param(object)))
setMethod("scale", "Cauchy", 
           function(x, center = TRUE, scale = TRUE) scale(param(x)))
           ### odd arg-list due to existing function in base package 

## wrapped replace methods
setMethod("location<-", "Cauchy", 
           function(object, value) 
               new("Cauchy", location = value, scale = scale(object)))
setMethod("scale<-", "Cauchy", 
           function(object, value) 
               new("Cauchy", location = location(object), scale = value))

## convolution operator for normal distributions

setMethod("+", c("Cauchy","Cauchy"),
          function(e1,e2){
              new("Cauchy", scale = scale(e1) + scale(e2), 
                   location = location(e1) + location(e2),  .withArith = TRUE)
          })

## some exact arithmetic methods
setMethod("+", c("Cauchy","numeric"),
          function(e1, e2){
              if (length(e2)>1) 
                  stop("length of operator must be 1")
              Cauchy(location = location(e1) + e2, scale = scale(e1))  
          })

setMethod("*", c("Cauchy","numeric"),
          function(e1, e2){
            if (length(e2)>1) 
                stop("length of operator must be 1")
            if(isTRUE(all.equal(e2,0))) 
               Dirac(location = 0) 
            else                
               Cauchy(location = location(e1) * e2, 
                      scale = scale(e1) * abs(e2))  
          })

