################################
##
## Class: BinomParameter
##
################################


## Access Methods
setMethod("size", "BinomParameter", function(object) object@size)
setMethod("prob", "BinomParameter", function(object) object@prob)
## Replace Methods
setReplaceMethod("size", "BinomParameter", 
                  function(object, value){ object@size <- value; object})
setReplaceMethod("prob", "BinomParameter", 
                  function(object, value){ object@prob <- value; object})


setValidity("BinomParameter", function(object){
  if(length(prob(object)) != 1)
    stop("prob has to be a numeric of length 1")    
  if(prob(object) < 0)
    stop("prob has to be in [0,1]")
  if(prob(object) > 1)
    stop("prob has to be in [0,1]")
  if(length(size(object)) != 1)
    stop("size has to be a numeric of length 1")    
  if(size(object) < 1)
    stop("size has to be a natural greater than 0")
  if(!identical(floor(size(object)), size(object)))
    stop("size has to be a natural greater than 0")    
  else return(TRUE)
})


################################
##
## Class: binomial distribution
##
################################

Binom <- function(size = 1,prob = 0.5) new("Binom", size = size, prob = prob)

## wrapped access methods
setMethod("prob", "Binom", function(object) prob(param(object)))
setMethod("size", "Binom", function(object) size(param(object)))
## wrapped replace methods
setMethod("prob<-", "Binom", 
           function(object, value) new("Binom", prob = value, 
                                        size = size(object)))
setMethod("size<-", "Binom", 
           function(object, value) new("Binom", prob = prob(object), 
                                        size = value))

## Convolution for two binomial distributions Bin(n1,p1) and Bin(n2,p2)
## Distinguish cases 
## p1 == p2 und p1 != p2


setMethod("+", c("Binom","Binom"),
          function(e1,e2){
            newsize <- size(e1) + size(e2)
            
            if(isTRUE(all.equal(prob(e1),prob(e2))))    
               return(new("Binom", prob = prob(e1), size = newsize, 
                          .withArith = TRUE))
            
            return(as(e1, "LatticeDistribution") + e2)
          })
