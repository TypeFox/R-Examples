"simulatesums" <-function(nrep, n, DIST, ...){
	temp <- matrix(DIST(n * nrep, ...), ncol = n, nrow = nrep)
	c(temp %*% rep(1., n))
	}