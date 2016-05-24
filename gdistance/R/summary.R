# Author: Jacob van Etten
# Date: September 2010
# Version 1.0
# Licence GPL v3


if (!isGeneric("summary")) {
	setGeneric("summary", function(object, ...)
		standardGeneric("summary"))
}	

setMethod('summary', signature(object='TransitionLayer'), 
	function(object, ...) {
		summary(transitionMatrix(object))
	}
)

setMethod('summary', signature(object='TransitionStack'), 
	function(object, ...) {
		n <- nlayers(object)
		result <- vector("list", length=n)
		for(i in 1:n)
		{
			result[i] <- summary(object@transition[i]@transitionMatrix)
		}
		return(result)
	}
)
