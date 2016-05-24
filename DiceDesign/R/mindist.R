mindist <- function(design){
	# Compute the Mindist criterion
	# A higher MINDIST corresponds to a more regular
	# scaterring of design points
	# J.FRANCO, 2006.10.23
	# input : design of n experiments

	X <- as.matrix(design)
	n <- dim(X)[1]
	# To check the experimental region
	if ( min(X)<0 || max(X)>1 ){
		warning("The design is rescaling into the unit cube [0,1]^d.")
		M <- apply(X,2,max)
		m <- apply(X,2,min)
		for (j in 1:dim(X)[2]){
			X[,j] <- (X[,j]-m[j])/(M[j]-m[j])
		}	
	}
	
	Distance <- as.matrix(dist(X))
	diag(Distance)<-1.0E30
	# Compute the min on the columns
	min1 <- apply(Distance,2,min)

	return(min(min1))
}

