`sumFunc` <-
function(x, cutoff) {
	tmp <- x/sum(x)
	final <- rep(0, length(x))
	for(i in 1:length(x)) {
		final[i] <- sum(tmp[1:i])
	}
	dim <- min(which(final > cutoff))
	return(dim)
}

