eds1 <- function(x, q=1, retq = TRUE) {
	bsums <- function(x) {
		if (all(x==0)) return(0)
		sum(-x*log(x),na.rm=TRUE)
	}
	x[x!=0] <- x[x!=0]^q
	if (!is.null(dim(x))) a <- apply(x,1,eds1,retq=retq)
	else {
		x <- matrix(x,nrow=1)
		rs <- sum(x)
		if (rs==0) stop("Zero data")
		x <- x/rs
		a <- bsums(x)
		if (retq) {
			a <- exp(a)
		}
	}
	a
}

