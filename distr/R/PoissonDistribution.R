
################################
##
## Class: PoisParameter
##
################################

## Access Methods
setMethod("lambda", "PoisParameter", function(object) object@lambda)
## Replace Methods
setReplaceMethod("lambda", "PoisParameter", 
                  function(object, value){ object@lambda <- value; object})

setValidity("PoisParameter", function(object){
  if(length(lambda(object)) != 1)
    stop("lambda has to be a numeric of length 1")    
  if(lambda(object) < 0)
    stop("lambda has to be not negative")
  else return(TRUE)
})

################################
##
## Class: Poisson distribution 
##
################################

Pois <- function(lambda = 1) new("Pois", lambda = lambda)

## wrapped access methods
setMethod("lambda", "Pois", function(object) lambda(param(object)))
## wrapped replace methods
setMethod("lambda<-", "Pois", 
           function(object, value) new("Pois", lambda = value))

setMethod("+", c("Pois","Pois"),
          function(e1,e2){
            new("Pois", lambda = lambda(e1) + lambda(e2), .withArith = TRUE)
          })
