`setParams` <-
function(x, s) {
	
	if (!is.gp(x) && !is.gp.list(x)) {
		stop("x must be of type gp or gp.list")
	} 
	if (length(s) != x$numDim){
		stop("s must be same length as x$numDim")
	}

	x$params = s	
	if (is.gp.list(x)) {
		for (i in 1:x$numGPs) {
			x[[1]]$params = s
		}
	}
	return (x)
}

