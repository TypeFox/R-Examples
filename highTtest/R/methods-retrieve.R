if(!isGeneric("CK")){
  setGeneric(name="CK", 
             def=function(x,...){standardGeneric("CK")})
}

setMethod("CK",
  signature = c(x="highTtest"),
  definition = function(x, ...){
    return(x@CK)
  }
)

if(!isGeneric("pi_alt")){
  setGeneric(name="pi_alt", 
             def=function(x,...){standardGeneric("pi_alt")})
}

setMethod("pi_alt",
  signature = c(x="highTtest"),
  definition = function(x, ...){
    return(x@pi1)
  }
)

if(!isGeneric("pvalue")){
  setGeneric(name="pvalue", 
             def=function(x,...){standardGeneric("pvalue")})
}

setMethod("pvalue",
  signature = c(x="highTtest"),
  definition = function(x, ...){
    return(x@pvalue)
  }
)

if(!isGeneric("ST")){
  setGeneric(name="ST", 
             def=function(x,...){standardGeneric("ST")})
}

setMethod("ST",
  signature = c(x="highTtest"),
  definition = function(x, ...){
    if(is.null(x@ST)){
      cat("ST method not calculated for provided object.\n")
      return(NULL)
    } else {
      return(x@ST)
    }
  }
)


if(!isGeneric("BH")){
  setGeneric(name="BH", 
             def=function(x,...){standardGeneric("BH")})
}

setMethod("BH",
  signature = c(x="highTtest"),
  definition = function(x, ...){
    if(is.null(x@BH)){
      cat("BH method not calculated for provided object.\n")
      return(NULL)
    } else {
      return(x@BH)
    }
  }
)

