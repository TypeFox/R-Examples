############################################################
# generating functions for BiasTypes
############################################################

symmetricBias <- function(name = "symmetric Bias")
     new("symmetricBias", name = name)

positiveBias <- function(name = "positive Bias")
     new("onesidedBias", name = name, sign = 1)

negativeBias <- function(name = "negative Bias")
     new("onesidedBias", name = name, sign = -1)
     
asymmetricBias <- function(name = "asymmetric Bias", nu = c(1,1))
     new("asymmetricBias", name = name, nu = nu)     
     
### Accessors & Replacements

setMethod("name", "BiasType", function(object) object@name)
setReplaceMethod("name", "BiasType", 
                  function(object, value) {object@name <- value; object})


setMethod("sign", "onesidedBias", function(x) x@sign)
setReplaceMethod("sign", "onesidedBias", 
                  function(object, value) 
                     {if (abs(trunc(value))!=1) stop("Left value has to be +-1")
                      object@sign <- value; object}
                  )

setMethod("nu", "asymmetricBias", function(object) object@nu)
setReplaceMethod("nu", "asymmetricBias", 
                  function(object, value) 
                     { if(!length(value)==2 || any(value>1) || 
                            any (value<=0) || max(value)<1 )   
                           stop("Left value has to be in (0,1]x(0,1] with maximum 1")
                       object@nu <- value; object}
                  )



