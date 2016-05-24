################################
##
## Class: UnifParameter
##
################################


## Access Methods
setMethod("Min", "UnifParameter", function(object) object@Min)
setMethod("Max", "UnifParameter", function(object) object@Max)
## Replace Methods
setReplaceMethod("Min", "UnifParameter", 
                  function(object, value){ object@Min <- value; object})
setReplaceMethod("Max", "UnifParameter", 
                  function(object, value){ object@Max <- value; object})


setValidity("UnifParameter", function(object){
  if(length(Min(object)) != 1)
    stop("Min has to be a numeric of length 1")    
  if(length(Max(object)) != 1)
    stop("Max has to be a numeric of length 1")    
  if(Min(object) >= Max(object))
    stop("Min has to be less than Max")
  else return(TRUE)
})


################################
##
## Class: uniform distribution
##
################################

Unif <- function(Min = 0, Max = 1) new("Unif", Min = Min, Max = Max)

## wrapped access methods
setMethod("Min", "Unif", function(object) Min(param(object)))
setMethod("Max", "Unif", function(object) Max(param(object)))

## wrapped replace methods
setMethod("Min<-", "Unif", 
           function(object, value) new("Unif", Min = value, Max = Max(object)))
setMethod("Max<-", "Unif", 
           function(object, value) new("Unif", Min = Min(object), Max = value))

## extra methods for Unif distribution
setMethod("+", c("Unif","numeric"),
          function(e1, e2){
            if (length(e2)>1) stop("length of operator must be 1")
            new("Unif", Min = Min(e1) + e2, Max = Max(e1) + e2, 
                .withArith = TRUE) 
          })

setMethod("*", c("Unif","numeric"),
          function(e1, e2){
            if (length(e2)>1) stop("length of operator must be 1")
            if(e2 == 0) return(new("Dirac", location = 0, .withArith = TRUE))
            if(e2 > 0) 
              new("Unif", Min = Min(e1) * e2, Max = Max(e1) * e2, 
                   .withArith = TRUE)
            else
              new("Unif", Min = Max(e1) * e2, Max = Min(e1) * e2, 
                  .withArith = TRUE)
          })

