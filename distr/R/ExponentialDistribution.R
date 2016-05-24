################################
##
## Class: ExpParameter
##
################################

## Access methods
setMethod("rate", "ExpParameter", function(object) object@rate)
## Replace methods
setReplaceMethod("rate", "ExpParameter", 
        function(object, value){ object@rate <- value; object})


setValidity("ExpParameter", function(object){
  if(length(rate(object)) != 1)
    stop("rate has to be a numeric of length 1")    
  if(rate(object) <= 0)
    stop("rate has to be positive")
  else return(TRUE) })


################################
##
## Class: exponential distribution
##
################################

Exp <- function(rate = 1) new("Exp", rate = rate)

## wrapped access methods
setMethod("rate", "Exp", function(object) rate(param(object)))

## wrapped replace methods
setMethod("rate<-", "Exp", function(object, value) new("Exp", rate = value))

setMethod("*", c("Exp","numeric"),
       function(e1, e2){
         if (length(e2)>1) stop("length of operator must be 1")
         if(isTRUE(all.equal(e2,0)))  
            return(new("Dirac", location = 0, .withArith = TRUE))
         if(e2 > 0) return(new("Exp", rate = rate(e1) / e2, .withArith = TRUE))
         return(-1 * as(Exp(rate = rate(e1) / abs(e2)), "AbscontDistribution"))
       })

setMethod("shape", "Exp", function(object)1)

