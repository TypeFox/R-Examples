# Expected values:
KE <- function(LAD, N=100, kneighbors=5, returnmatrix=FALSE){

	bound <- (N/LAD)^(1/3)
	
	x <- runif(N,0,bound)
	y <- runif(N,0,bound)
	z <- runif(N,0,bound)
	nleaf <- N
	res <- matrix(nrow=N,ncol=N)
	res[] <- 0
	
	f <- .Fortran("Kn", as.double(x),as.double(y),as.double(z),as.integer(nleaf),as.double(res),
		PACKAGE="YplantQMC")
	
	m <- f[[5]]
	m <- matrix(m, ncol=nleaf)
	
	if(returnmatrix)return(m)
	
	m[m == 0] <- NA
	meanminn <- function(x,k)mean(sort(x)[1:k], na.rm=TRUE)
	
	if(length(kneighbors) == 1){
		o <- apply(m, 1, meanminn,kneighbors)
		return(mean(o))
	} else {
		o <- list()
		for(j in 1:length(kneighbors)){
			o[[j]] <- apply(m, 1, meanminn,kneighbors[j])
		}
		return(sapply(o,mean))
	}
return(mean(o))
}

