if(!isGeneric("deviance"))
  setGeneric(name = "deviance", def = function(object, ...) standardGeneric("deviance"))

## method deviance form models of class glimML (betabin and negbin)
setMethod(f = "deviance", signature = "glimML", definition =  function(object, ...) object@dev)
