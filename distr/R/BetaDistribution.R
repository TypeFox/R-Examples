
################################
##
## Class: BetaParameter
##
################################


## Access methods

setMethod("shape1", "BetaParameter", function(object) object@shape1)
setMethod("shape2", "BetaParameter", function(object) object@shape2)
setMethod("ncp", "BetaParameter", function(object) object@ncp)

## Replace Methoden

setReplaceMethod("shape1", "BetaParameter", 
                  function(object, value){ object@shape1 <- value; object})
setReplaceMethod("shape2", "BetaParameter", 
                  function(object, value){ object@shape2 <- value; object})
setReplaceMethod("ncp", "BetaParameter", 
                  function(object, value){ object@ncp <- value; object})


setValidity("BetaParameter", function(object){
  if(length(shape1(object)) != 1)
    stop("shape1 has to be a numeric of length 1")    
  if(shape1(object) <= 0)
    stop("shape1 has to be positive")
  if(length(shape2(object)) != 1)
    stop("shape2 has to be a numeric of length 1")    
  if(shape2(object) <= 0)
    stop("shape2 has to be positive")
  if(length(ncp(object)) != 1)
    stop("ncp has to be a numeric of length 1")      
  else return(TRUE)
}
)


################################
##
## Class: BetaDistribution
##
################################

Beta <- function(shape1 = 1, shape2 = 2, ncp = 0) 
                 new("Beta", shape1 = shape1, shape2 = shape2, ncp = ncp)

## wrapped access methods
setMethod("shape1", "Beta", function(object) shape1(param(object)))
setMethod("shape2", "Beta", function(object) shape2(param(object)))
setMethod("ncp", "Beta", function(object) ncp(param(object)))

## wrapped replace methods
setMethod("shape1<-", "Beta", 
           function(object, value) new("Beta", shape1 = value, 
                                    shape2 = shape2(object), ncp = ncp(object)))
setMethod("shape2<-", "Beta", 
           function(object, value) new("Beta", shape1 = shape1(object), 
                                        shape2 = value, ncp = ncp(object)))
setMethod("ncp<-", "Beta", 
           function(object, value) new("Beta", shape1 = shape1(object), 
                                        shape2 = shape2(object), ncp = value))

setMethod("-", c("numeric","Beta"),
           function(e1, e2) {if(isTRUE(all.equal(e1,1))&&
                                isTRUE(all.equal(ncp(e2),0)))
                             return(Beta(shape1=shape2(e2),shape2=shape1(e2)))
                             else e1-as(e2,"AbscontDistribution")})
