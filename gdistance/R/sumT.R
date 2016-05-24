sumReciprocal <- function(x1, x2)
{
	if(matrixValues(x1) == "conductance") 
	{
		x1@transitionMatrix@x <- 1 / x1@transitionMatrix@x
		x1@transitionMatrix@x[x1@transitionMatrix@x == Inf] <- 0
	}
	if(matrixValues(x2) == "conductance") 
	{
		x2@transitionMatrix@x <- 1 / x2@transitionMatrix@x
		x2@transitionMatrix@x[x2@transitionMatrix@x == Inf] <- 0
	}
	newTransition <- x1 + x2
	newTransition@transitionMatrix@x <- 1 / newTransition@transitionMatrix@x
	if(sum(newTransition@transitionMatrix@x == Inf) > 0) 
	{
		warning("Inf values introduced. Set to 0.")
		newTransition@transitionMatrix@x[newTransition@transitionMatrix@x == Inf] <- 0
	}
	matrixValues(newTransition) <- "conductance"
	return(newTransition)
}
