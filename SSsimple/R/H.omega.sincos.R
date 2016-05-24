H.omega.sincos <-
function(x, u) {
	n <- length(x)
	d <- length(u)
	y <- rep(x, d)
	w <- rep(u, each=n)
	H <- matrix( NA, n, d*2 )
	H[ , seq(1, d*2, by=2)] <- sin(y*w)
	H[ , seq(2, d*2, by=2)] <- cos(y*w)
	return(H)
}

