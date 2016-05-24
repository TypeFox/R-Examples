coverage <- function(design){
	# Compute the coverage measure
	# For a regular mesh, c=0
	# input : design of n experiments

	X <- as.matrix(design)
	n <- dim(X)[1]
	dimension <- dim(X)[2]

	if ( n < dimension ){
    		stop('Warning : the number of points is lower than the dimension.')
	}
	# To check the experimental region
	if ( min(X)<0 || max(X)>1 ){
		warning("The design is rescaling into the unit cube [0,1]^d.")
		M <- apply(X,2,max)
		m <- apply(X,2,min)
		for (j in 1:dim(X)[2]){
			X[,j] <- (X[,j]-m[j])/(M[j]-m[j])
		}	
	}
	
	# Compute the distance between the points.
	Distance <- as.matrix(dist(X))
	diag(Distance) <- 10e3
	Dmin <- as.matrix(apply(Distance,2,min))
	gammabar <-(1/n)*sum(Dmin)
	s <- 0

	for (i in 1:n){
    		s <- s+ (Dmin[i]-gammabar)*(Dmin[i]-gammabar)
	}
	cov <-(1/gammabar)*((1/n)*s)^(1/2)
	return(cov)
}