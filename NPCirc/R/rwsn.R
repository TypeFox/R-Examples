rwsn<-function(n,xi,eta,lambda){
	if (!is.numeric(n)) stop("argument 'n' must be numeric")
	if (missing(xi) || length(xi)!=1 || !is.numeric(xi))
    	stop("The mean direction parameter 'xi' is mandatory. It must be numeric and have length 1")
 	if (missing(eta) || length(eta)!=1 || !is.numeric(eta))
    	stop("The scale parameter 'eta' is mandatory. It must be numeric and have length 1")
  	if (missing(lambda) || length(lambda)!=1 || !is.numeric(lambda))
    	stop("The shape parameter 'lambda' is mandatory. It must be numeric and have length 1")
    	xi <- conversion.circular(xi, units="radians", zero=0, rotation="counter")
	attr(xi, "class") <- attr(xi, "circularp") <- NULL
 	u1 <- rnorm(n)
    	u2 <- rnorm(n)
    	id <- (u2 > lambda * u1)
    	u1[id] <- (-u1[id])
    	y <- xi + eta * u1
	x <- y%%(2*pi)
	x <- circular(x)
	return(x)
}


