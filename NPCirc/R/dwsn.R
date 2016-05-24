dwsn<- function(x,xi,eta,lambda,K=NULL,min.k=20){
	x <- conversion.circular(x, units="radians", zero=0, rotation="counter")
  	xi <- conversion.circular(xi, units="radians", zero=0, rotation="counter")
  	attr(x, "class") <- attr(x, "circularp") <- NULL
  	attr(xi, "class") <- attr(xi, "circularp") <- NULL


	if (missing(xi) || length(xi)!=1 || !is.numeric(xi))
    	stop("The mean direction parameter 'xi' is mandatory. It must be numeric and have length 1")
 	if (missing(eta) || length(eta)!=1 || !is.numeric(eta))
    	stop("The scale parameter 'eta' is mandatory. It must be numeric and have length 1")
  	if (missing(lambda) || length(lambda)!=1 || !is.numeric(lambda))
    	stop("The shape parameter 'lambda' is mandatory. It must be numeric and have length 1")
	x <- x[!is.na(x)]
	n <- length(x)
	if (sum(is.na(x))>0) warning("Missing values were removed")
	if (n==0) stop("No observations (at least after removing missing values)")
  	if (is.null(K)){
    		range <- abs(xi-x)
    		K <- (range+6*eta)%/%(2*pi)+1
    		K <- max(min.k, K)
  	}else{
		if (!is.numeric(K) | K<=0){
			warning("Argument 'K' must be a positive integer. 'K=min.k' was used")	
			K<-min.k
		}
	}
	fx<-numeric(n)
	for (i in 1:n){
		val <- (x[i]+2*pi*seq(-K,K,1)-xi)/eta
		suma <- sum(dnorm(val)*pnorm(lambda*val))
		fx[i] <- 2/eta*suma
	}
	return(fx)
}


