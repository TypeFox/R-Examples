.cfilter <- function(x, alpha=0.20, side="both"){
	n <- length(x)
	x.order <- order(x)  
	x <- sort(x)
  
	## upper
	q.up <- quantile(x, prob=1-alpha)
	x.up <- x[which( x > q.up)] - q.up
	n.up <- length(x.up)
	rate.up <- log(2)/median(x.up)
	dn.up <- 0
	if( sum(x.up > 1/rate.up) > 0 )
		dn.up <- max(pmax( pexp( x.up, rate.up) - (1:n.up - 1)/n.up, 0)[ which( x.up > 1/rate.up)])
	n0.up <- round(dn.up*n.up)

	## lower
	q.low <- quantile(x, prob=alpha)
	x.low <- rev(q.low - x[which( x < q.low)])
	n.low <- length(x.low)
	rate.low <- log(2)/median(x.low)
	dn.low <- 0
	if( sum(x.low > 1/rate.low) > 0 )
		dn.low <- max(pmax( pexp( x.low, rate.low) - (1:n.low - 1)/n.low, 0)[ which( x.low > 1/rate.low)])
	n0.low <- round(dn.low*n.low)

  	## filter
	x.na <- switch(side, 
		both={if(n0.up > 0) x[ (n - n0.up + 1):n ] <- NA;
			if(n0.low > 0) x[ 1:n0.low ] <- NA;
			x},
		upper={if(n0.up > 0) x[ (n - n0.up + 1):n ] <- NA;
			x},
		lower={if(n0.low > 0) x[ 1:n0.low ] <- NA;
			x})
	x.na <- x.na[order(x.order)]
	return( x.na )
}

.cfilter.pert <- function(x, m=5, gamma=0.5, cutoff=0.5, alpha=0.20, side="both", seed){
	if( missing(seed) ) seed <- 1000
	if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
		seed.keep <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
		on.exit(assign(".Random.seed", seed.keep, envir = .GlobalEnv))
	}
	set.seed(seed)
	n <- length(x)
	x.na <- matrix(NA, m, n)
	s <- mad(x)
	x.na[1,] <- .cfilter(x, alpha, side)
	for(i in 2:m){
		xpert <- x + gamma*s*rnorm(n)
		x.na[i,] <- .cfilter(xpert, alpha, side)
	}
	x[ colMeans(is.na(x.na)) > cutoff] <- NA
	return(x)
}

.cfilter.iter <- function(x, alpha=0.20, miter=3){
	iter <- 0
	n <- length(x)
	id <- 1:n
	x.old <- x
	converge <- 0
	while( converge == 0 & iter < miter ){
		iter <- iter + 1
		x <- .cfilter.pert( x, alpha=alpha)
		id <- id[ !is.na(x) ] 
		if( !any(is.na(x)) ) converge <- 1
		x <- na.omit(x)
	}
	x.na <- rep(NA, n)
	x.na[id] <- x.old[id]
	return( x.na)
}



