dlogr.predict <-
function (x,BetaHat) {
# function to get predicted values for a set of 1 or more x's
	if (!is.matrix(x)) {
	   dim(x) <- c(1,length(x))
	}
	return(1/(1+exp(-x %*% BetaHat)))
}

