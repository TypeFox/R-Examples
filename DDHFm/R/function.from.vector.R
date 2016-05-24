"function.from.vector" <-
function(x, y, argument.vect) {

	#	x must be in ascending order
	#	the lengths of x and y must match
	#	argument.vect can be a vector
	
	indices <- sapply(argument.vect, which.min.diff, x)
	return(y[indices])
}

