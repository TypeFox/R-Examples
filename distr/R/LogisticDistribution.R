
################################
##
## Class: LogisParameter
##
################################
## Access Methods
setMethod("location", "LogisParameter", function(object) object@location)
setMethod("scale", "LogisParameter", 
           function(x, center = TRUE, scale = TRUE) x@scale)
## Replace Methods
setReplaceMethod("location", "LogisParameter", 
                  function(object, value){ object@location <- value; object})
setReplaceMethod("scale", "LogisParameter", 
                  function(object, value){ object@scale <- value; object})

setValidity("LogisParameter", function(object){
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
## Class: logistic distribution
##
################################

Logis <- function(location = 0, scale = 1) 
         new("Logis", location = location, scale = scale)

## wrapped access methods
setMethod("location", "Logis", function(object) location(param(object)))
setMethod("scale", "Logis", 
           function(x, center = TRUE, scale = TRUE) scale(param(x)))

## wrapped replace methods
setMethod("location<-", "Logis", function(object, value) 
           new("Logis", location = value, scale = scale(object)))
setMethod("scale<-", "Logis", function(object, value) 
           new("Logis", location = location(object), scale = value))

## extra Methoden for Logistic distribution
setMethod("+", c("Logis","numeric"),
          function(e1, e2){
            if (length(e2)>1) stop("length of operator must be 1")
            new("Logis", location = location(e1) + e2, 
                 scale = scale(e1), .withArith = TRUE) 
          })

setMethod("*", c("Logis","numeric"),
          function(e1, e2){
            if (length(e2)>1) stop("length of operator must be 1")
            if(isTRUE(all.equal(e2,0)))  
               return(new("Dirac", location = 0, .withArith = TRUE))
            new("Logis", location = location(e1) * e2, 
                 scale = scale(e1) * abs(e2), .withArith = TRUE)
          })


