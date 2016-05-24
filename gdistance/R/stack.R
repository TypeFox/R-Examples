# Author: Jacob van Etten jacobvanetten@yahoo.com
# IE University
# Date :  June 2010
# Version 1.0
# Licence GPL v3

setMethod("stack", signature(x="TransitionLayer"), function(x, ...) 
	{
		newStack <- as(x, "TransitionStack")
		objectList <- list(...)
		TData <- .createTData(x, objectList)
		newStack@transition <- TData
		return(newStack)	
	} 
)

setMethod("stack", signature(x='TransitionStack'), function(x, ...) 
	{
		newStack <- as(x, "TransitionStack")
		objectList <- list(...)
		TData <- .createTData(x, objectList)
		newStack@transition <- TData
		return(newStack)	
	} 
)

.createTData <- function(x, objectList)
{
	nobj <- length(objectList)
	if(nobj<1) {stop("more than one object is needed to stack")}
	TD <- transitionData(x)
	for(i in 1:nobj) {TD <- c(TD,transitionData(objectList[[i]]))}
	return(TD)
}