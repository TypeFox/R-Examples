`setGPNames` <-
function(x, s) {
	if (!is.gp.list(x)) { 
		stop("error: x must be of type gp.list")
	}
	if (length(s) != x$numGPs) {
		stop("s must be same length as x$numGPs")
	}
	x$names = s
	return (x)
}

