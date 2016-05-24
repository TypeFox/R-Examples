################################
##
## Class: HyperParameter
##
################################


## Access Methods
setMethod("m", "HyperParameter", function(object) object@m)
setMethod("n", "HyperParameter", function(object) object@n)
setMethod("k", "HyperParameter", function(object) object@k)
## Replace Methods
setReplaceMethod("m", "HyperParameter", 
                  function(object, value){ object@m <- value; object})
setReplaceMethod("n", "HyperParameter", 
                  function(object, value){ object@n <- value; object})
setReplaceMethod("k", "HyperParameter", 
                  function(object, value){ object@k <- value; object})




setValidity("HyperParameter", function(object){
  if(length(m(object)) != 1)
    stop("m has to be a numeric of length 1")    
  if(m(object) < 0)
    stop("m has to be a not negative natural")
  if(length(n(object)) != 1)
    stop("n has to be a numeric of length 1")    
  if(n(object) < 0)
    stop("n has to be a not negative natural")
  if(length(k(object)) != 1)
    stop("k has to be a numeric of length 1")    
  if(k(object) < 0)
    stop("k has to be a not negative natural")
  
  if(!identical(floor(m(object)), m(object)))
    stop("m has to be a not negative natural")
  if(!identical(floor(n(object)), n(object)))
    stop("n has to be a not negative natural")
  if(!identical(floor(k(object)), k(object)))
    stop("k has to be a not negative natural")
  if(k(object) > m(object) + n(object))
    stop("k has to be less or equal than m + n")    
  else return(TRUE)
})



################################
##
## Class: hypergeometric distribution
##
################################

Hyper <- function(m = 1, n = 1, k = 1) new("Hyper", m = m, n = n, k = k)

## wrapped access methods
setMethod("m", "Hyper", function(object) m(param(object)))
setMethod("n", "Hyper", function(object) n(param(object)))
setMethod("k", "Hyper", function(object) k(param(object)))
## wrapped replace methods
setMethod("m<-", "Hyper", 
           function(object, value) 
                    new("Hyper", m = value, n = n(object), k = k(object)))
setMethod("n<-", "Hyper", 
           function(object, value) 
                    new("Hyper", m = m(object), n = value, k = k(object)))
setMethod("k<-", "Hyper", 
           function(object, value) 
                    new("Hyper", m = m(object), n = n(object), k = value))
