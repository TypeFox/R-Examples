"simulatereg" <-function(x, nrep, alpha = 0., beta = 1., DIST = rnorm, ...){
	n <- length(x)
	X <- matrix(x, ncol = n, nrow = nrep, byrow = TRUE)
	Y <- (X * beta + alpha) + matrix(DIST(nrep * n, ...), ncol = n, nrow = nrep)
	xc <- x - mean(x)
	ssxx <- sum(xc^2.)
	ssxy <- (Y * rep(xc, rep(nrep, n))) %*% rep(1., n)
	my <- Y %*% rep(1./n, n)
	ssyy <- ((Y - rep(my, n))^2.) %*% rep(1., n)
	b <- c(ssxy/ssxx)
	a <- c(my - b * mean(x))
	s2 <- c((ssyy - b * ssxy)/(n - 2.))
	#list(my = my, Y = Y, a = b, b = b, s2 = s2, ssxx = ssxx, ssyy = ssyy, 
	#ssxy = ssxy)
	list(a = a, b = b, s = sqrt(s2), ssx = ssxx, df = n - 2., x = x)
	}