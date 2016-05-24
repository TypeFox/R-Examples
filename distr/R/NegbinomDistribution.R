################################
##
## Class: NbinomParameter
##
################################

## Access Methods
setMethod("size", "NbinomParameter", function(object) object@size)
setMethod("prob", "NbinomParameter", function(object) object@prob)
## Replace Methods
setReplaceMethod("size", "NbinomParameter", 
                  function(object, value){ object@size <- value; object})
setReplaceMethod("prob", "NbinomParameter", 
                  function(object, value){ object@prob <- value; object})

setValidity("NbinomParameter", function(object){
  if(length(prob(object)) != 1)
    stop("prob has to be a numeric of length 1")    
  if(prob(object) <= 0)
    stop("prob has to be in (0,1)")
  if(prob(object) >= 1)
    stop("prob has to be in (0,1)")
  if(length(size(object)) != 1)
    stop("size has to be a numeric of length 1")    
  if(size(object) < 0)
    stop("size has to be non-negative")
#  if(!identical(floor(size(object)), size(object)))
#    stop("size has to be a not negative natural")
  else return(TRUE)
})


################################
##
## Class: negative binomial distribution
##
################################

Nbinom <- function(size = 1,prob = 0.5) new("Nbinom", size = size, prob = prob)

## wrapped access methods
setMethod("prob", "Nbinom", function(object) prob(param(object)))
setMethod("size", "Nbinom", function(object) size(param(object)))
## wrapped replace methods
setMethod("prob<-", "Nbinom", 
           function(object, value) 
               new("Nbinom", prob = value, size = size(object)))
setMethod("size<-", "Nbinom", 
           function(object, value) 
               new("Nbinom", prob = prob(object), size = value))

setMethod("+", c("Nbinom","Nbinom"),
          function(e1,e2){
            newsize <- size(e1) + size(e2)
            
            if(isTRUE(all.equal(prob(e1),prob(e2))))    
              return(new("Nbinom", size = newsize, prob = prob(e1), 
                         .withArith = TRUE))           
            return(as(e1, "LatticeDistribution") + e2)
          })
