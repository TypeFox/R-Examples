.gy.filt <- function(x, alpha, it=TRUE){
	x <- as.matrix(x)
	xs <- scale(x, apply(x, 2, median, na.rm=TRUE), apply(x, 2, mad, na.rm=TRUE))
	xs2 <- xs^2
	xs2.na <- switch( it + 1, 
		apply(xs2, 2, .gy.filt.uni, alpha=alpha),
		apply(xs2, 2, .gy.filt.uni.it, alpha=alpha))
	x.na <- x
	x.na[ which(is.na(xs2.na)) ] <- NA
	if( ncol(x.na) == 1) x.na <- c(x.na)
	return( x.na )
}


.gy.filt.uni <- function(v, alpha){
	n <- length(v)
	v.order <- order(v)	
	v <- sort(v)
	i0 <- which(v < qchisq( alpha, 1 ))
	n0 <- 0
	if( length(i0) > 0){
		i0 <- rev(i0)[1]
		dn <- max( pmax( pchisq( v[i0:n], 1) - (i0:n - 1)/n, 0)) 
		n0 <- round(dn*n)
    		## newly added: 2015-04-09
		##in0 <- i0 + which.max(pmax(pchisq(v[i0:n], 1) - (i0:n - 1)/n, 0)) - 1
	} 
	v <- v[ order(v.order) ] 
	v.na <- v
	if(n0 > 0) v.na[ v.order[ (n - n0 + 1):n] ] <- NA
	## newly added: 2015-04-09
  	#if(n0 > 0) v.na[ v.order[ in0:(in0 + n0 - 1) ] ] <- NA
	return( v.na )
}


.gy.filt.uni.it <- function(v, alpha, miter=10){
	converge <- 0	
	iter <- 0
	n <- length(v)
	id <- 1:n
	v.old <- v
	while( converge == 0 & iter < miter ){
		iter <- iter + 1
		v <- .gy.filt.uni( v, alpha )
		id <- id[ !is.na(v) ] 
		if( !any(is.na(v)) ) converge <- 1
		v <- na.omit(v)
	}
	v.out <- rep(NA, n)
	v.out[id] <- v
	return( v.out)
}

.impute.simple <- function(x, loc){
	u <- is.na(x)
	x[is.na(x)] <- 0
	x <- x + sweep(u, 2, loc, "*")
	x
}


TSGS <- function(x, alpha=0.95, it=TRUE, init=c("emve","sign","qc","huber","imputed"), partial.impute, ...){
	xcall <- match.call()
	init <- match.arg(init)
	
	## check dat
	if(is.data.frame(x) | is.matrix(x))
		x <- data.matrix(x)
	else stop("Data matrix must be of class matrix or data.frame.")
	if(any(is.na(x))) warning("Data matrix contains missing values.")
	## June 12, 2012
	## Only allow up to p=50
	n <- nrow(x)
	p <- ncol(x)
	if( p >200 | p < 2 ) stop("Column dimension of 'x' must be in between 2 and 200.")

	## filter step
	xf <- .gy.filt(x, alpha=alpha, it=it)

	## partial imputation step 
	xf_pi <- xf
	if( missing(partial.impute) ) partial.impute <- FALSE 
	if( partial.impute ){
		xmed <- apply(x, 2, median)
		ximp <- .impute.simple(xf_pi, xmed)
		aid <- which(rowSums(!is.na(xf_pi)) == p)
		uid <- which(rowSums(!is.na(xf_pi)) < p)
		n0 <- n/2 + (p+1)
		if( n0 > length(aid) ){
			fid <- sample( uid, n0 - length(aid))
			xf_pi[fid,] <- ximp[fid,] 
		}
	}
	
	res <- GSE(xf_pi, init=init, ...)
	res <- new("TSGS",
		call = xcall,
		S = res@S,
		mu = res@mu,
		xf = xf,
		sc = res@sc,
		mu0 = res@mu0,
		S0 = res@S0, 
		iter = res@iter,
		eps = res@eps,
		estimator = "2SGS", 
		x = x,
		ximp = res@ximp,
		weights = res@weights,
		weightsp = res@weightsp,
		pmd = res@pmd,
		pmd.adj = res@pmd.adj,
		p = res@p,
		pu = res@pu)
	res
}
