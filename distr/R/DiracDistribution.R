#################################
###
### Class: DiracParameter
###
#################################


### Access Methods
setMethod("location", "DiracParameter", function(object) object@location)
### Replace Methods
setReplaceMethod("location", "DiracParameter", 
                  function(object, value){ object@location <- value; object})

validDiracParameter <- function(object){
  if(length(location(object)) != 1)
    stop("location has to be a numeric of length 1")    
  else return(TRUE)
}

setValidity("DiracParameter", validDiracParameter)



#################################
###
### Class: Dirac distribution
###
#################################

Dirac <- function(location = 0) new("Dirac", location = location)

### wrapped access methods
setMethod("location", "Dirac", function(object) location(param(object)))
### wrapped replace methods
setMethod("location<-", "Dirac", 
           function(object, value) new("Dirac", location = value))


## some exact arithmetic methods
setMethod("+", c("Dirac","numeric"),
          function(e1, e2){
            if (length(e2)>1) stop("length of operator must be 1")
            new("Dirac", location = location(e1) + e2, .withArith = TRUE)  
          })

setMethod("*", c("Dirac","numeric"),
          function(e1, e2){
            if (length(e2)>1) stop("length of operator must be 1")
            new("Dirac", location = location(e1) * e2, .withArith = TRUE)  
          })

setMethod("+", c("Dirac","Dirac"),
          function(e1,e2){
            new("Dirac", location = location(e1) + location(e2), 
                 .withArith = TRUE)
          })

setMethod("*", c("Dirac","Dirac"),
          function(e1,e2){
            new("Dirac", location = location(e1) * location(e2), 
                 .withArith = TRUE)
          })

setMethod("-", c("Dirac","Dirac"),
          function(e1,e2){
            new("Dirac", location = location(e1) - location(e2), 
                 .withArith = TRUE)
          })

setMethod("/", c("Dirac","Dirac"),
          function(e1,e2){
            if(isTRUE(all.equal(location(e2),0))) 
                stop("location parameter of divisor must not be 0")
            new("Dirac", location = location(e1) / location(e2), 
                .withArith = TRUE)
          })

setMethod("/", c("numeric","Dirac"),
          function(e1,e2){
            if(isTRUE(all.equal(location(e2),0))) 
                stop("location parameter of divisor must not be 0")
            new("Dirac", location = e1 / location(e2), 
                .withArith = TRUE)
          })

setMethod("+", c("Dirac","UnivariateDistribution"),
          function(e1,e2){
            location(e1) + e2
          })

setMethod("*", c("Dirac","UnivariateDistribution"),
          function(e1,e2){
            location(e1) * e2
          })

setMethod("+", c("UnivariateDistribution","Dirac"),
          function(e1,e2){
            location(e2) + e1
          })

setMethod("*", c("UnivariateDistribution","Dirac"),
          function(e1,e2){
            location(e2) * e1
          })
