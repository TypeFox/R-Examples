
HuberPairwise <- function( x, psi=c("huber","sign"), c0=1.345, computePmd=TRUE){
	## argument checks
	## check choices of psi
	psi <- match.arg(psi)

	## check x
	if(is.data.frame(x) | is.matrix(x))
		x <- data.matrix(x)
	else stop("Data matrix must be of class matrix or data.frame.")

	xcall <- match.call()

	## drop all rows with missing values (!!) :
	x_nonmiss <- is.na(x)*-1 + 1
	pp <- rowSums(x_nonmiss)
	pp_col <- colSums(x_nonmiss)
	ok <- which(pp > 0); not.ok <- which(pp == 0)
	if( length(not.ok) > 0) cat("Observations (rows): ", paste(not.ok, collapse=", "), 
		"\nare completely missing and will be dropped out from the estimation.\n")
	x.orig <- x
	x <- x[ ok,]
	
	## Cannot contain all obs with completely missing rows!!
	if( all(pp == 0) ) stop("All observations have missing values!")
	if( any(pp_col == 0) )stop("Data matrix cannot contain column(s) with completely missing data!")

	## Check c0
	if( c0 == 0 ) psi <- "sign"
	if( c0 < 0 ) stop( "Input 'c0' must be positive")

	## initalization
	n <- nrow(x); p <- ncol(x)
	if( p <= 1 ) stop("Dimension of data must be >= 2!")

	if(n <= p + 1)
		stop(if (n <= p) "n <= p -- you can't be serious!" else "n == p+1  is too small sample size")
	if(n < 2 * p)
		## p+1 < n < 2p
		warning("n < 2 * p, i.e., possibly too small sample size")

	## Compute univariate location and scale
	sc <- apply(x, 2, mad, na.rm=T)
	if( any(sc == 0) ) stop("Some scale for certain variable(s) is 0!")
	loc <- apply(x, 2, median, na.rm=T)

	## replace missing data with the corresponding column wise median
	a <- is.na(x)%*%diag(loc)
	x[is.na(x)] <- 0
	x <- x + a
	
	## Transformation
	x.tmp <- sweep(sweep(x, 2, loc, "-"), 2, sc, "/")
	x.tmp <- switch( psi, 
		huber=apply( x.tmp, 2, function(x) pmin( pmax( x, -c0), c0) ),
		sign=apply( x.tmp, 2, sign ) )
	
	## Compute correlation matrix on transformed data
	R <- matrix(NA, p, p)
	x.tmp.sig <- sqrt( apply(x.tmp, 2, var, na.rm=T)*(n-1)/n )
	for(i in 1:(p-1))
		for(j in (i+1):p)
			R[i,j] <- R[j,i] <- mean( x.tmp[,i] * x.tmp[,j], na.rm=T )/(x.tmp.sig[i]*x.tmp.sig[j])	
	diag(R) <- 1

	## Currently, we don't adjust for bias here
	## please see QC investigation folder for the full R code

	## Form covariance matrix
	D <- diag( sc )
	covariance <- D %*% R %*% D
	cov.chol <- tryCatch( chol(covariance), error=function(e) NA)
	if( !is.matrix(cov.chol) )  warning("Estimated covariance matrix is not positive definite. May consider increase the sample size of the data.")		

	estimator <- switch( psi, 
			huber="Huber Pairwise",
			sign="Quadrant Covariance")
	
	## compute pmd
	pmd <- pmd.adj <- rep(NA, nrow(x.orig))
	pu <- rowSums( !is.na(x.orig) )
	if( computePmd ){
		pmd.obj <- partial.mahalanobis( x.orig, loc, covariance)
		pmd <- pmd.obj@pmd
		pmd.adj <- pmd.obj@pmd.adj
	}
	
	res <- new("HuberPairwise",
		call = xcall,
		S = covariance,
		R = R,
		mu = loc,
		estimator = estimator, 
		x = x.orig,
		pmd = pmd,
		pmd.adj = pmd.adj,
		p = p,
		pu = pu)
	res		
}

