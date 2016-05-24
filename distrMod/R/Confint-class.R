setMethod("confint", signature(object="Confint", method="missing"),
           function(object, method) object@confint)

setMethod("call.estimate", signature(object="Confint"),
           function(object) object@call.estimate)
setMethod("name.estimate", signature(object="Confint"),
           function(object) object@name.estimate)
setMethod("completecases.estimate", signature(object="Confint"),
           function(object) object@completecases.estimate)
setMethod("samplesize.estimate", signature(object="Confint"),
           function(object, onlycompletecases = TRUE)
  	    (object@samplesize.estimate+
  	     (1-onlycompletecases)*sum(object@completecases.estimate==FALSE)))
setMethod("nuisance.estimate", signature(object="Confint"),
           function(object) object@nuisance.estimate)
setMethod("trafo.estimate", signature(object="Confint"),
           function(object) object@trafo.estimate)
setMethod("fixed.estimate", signature(object="Confint"),
           function(object) object@fixed.estimate)
setMethod("type", signature(object="Confint"),
           function(object) object@type)



