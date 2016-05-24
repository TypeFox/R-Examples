## access methods
setMethod("loc", "GumbelParameter", function(object) object@loc)
setMethod("scale", "GumbelParameter", 
    function(x, center = TRUE, scale = TRUE) x@scale)

## replace Methods
setReplaceMethod("loc", "GumbelParameter", 
    function(object, value){ object@loc <- value; object })
setReplaceMethod("scale", "GumbelParameter", 
    function(object, value){ object@scale <- value; object})


## generating function
Gumbel <- function(loc = 0, scale = 1){ new("Gumbel", loc = loc, scale = scale) }

## wrapped access methods
setMethod("loc", "Gumbel", function(object) loc(object@param))
setMethod("scale", "Gumbel", 
    function(x, center = TRUE, scale = TRUE) scale(x@param))

## wrapped replace methods
setMethod("loc<-", "Gumbel", 
    function(object, value){ 
        new("Gumbel", loc = value, scale = scale(object))
    })
setMethod("scale<-", "Gumbel", 
    function(object, value){ 
        if(length(value) != 1 || value <= 0)
            stop("'value' has to be a single positive number")
        new("Gumbel", loc = loc(object), scale = value)
    })

## extra methods for Gumbel distribution
setMethod("+", c("Gumbel","numeric"),
          function(e1, e2){
            if (length(e2)>1) stop("length of operator must be 1")
            new("Gumbel", loc = loc(e1) + e2, scale = scale(e1)) 
          })

setMethod("*", c("Gumbel","numeric"),
          function(e1, e2){
            if (length(e2)>1) stop("length of operator must be 1")
            if (isTRUE(all.equal(e2,0))) 
                return(new("Dirac", location = 0, .withArith = TRUE))
            new("Gumbel", loc = loc(e1) * e2, scale = scale(e1)*abs(e2))
          })

### Euler Mascheroni constant:
EULERMASCHERONICONSTANT <- -digamma(1) ### after http://mathworld.wolfram.com/Euler-MascheroniConstant.html (48)

### Apéry constant
##local helper function:
.fctApery <- function(n) (-1)^n*choose(2*n,n)*n^3
##
APERYCONSTANT <- -sum(sapply(1:50,.fctApery)^(-1))*5/2 ## after http://mathworld.wolfram.com/AperysConstant.html (8)
