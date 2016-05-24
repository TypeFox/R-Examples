
tmvnsim <- function(nsamp, k, lower=rep(-Inf,k), upper=rep(Inf,k), imod=rep(FALSE, k), means=rep(0, k), sigma=diag(1, k))
{
	QR <- FALSE 
	if(QR)
	{
#		Doesn't seem to work.
		evd <- eigen(sigma)
		tol <- sqrt(.Machine$double.eps)
		Positive <- (evd$values > max(tol * evd$values[1], 0))
		npos <- sum(Positive)
		sroot1 <- evd$vectors[, 1:npos] %*% diag(sqrt(evd$values[1:npos]), npos) %*% t(evd$vectors[, 1:npos])

		sroot <- round(t(qr.R(qr(sroot1))), 6)
		ord <- 1:k
	} else 
	{
		Q <- try(chol(sigma, pivot=TRUE), silent=TRUE)
		if(inherits(Q, "try-error")) stop("Error in Cholesky Decomposition")
		sroot <- round(t(Q), 6)
		pivot <- attr(Q, "pivot")
		ord <- order(pivot)
		lower=lower[pivot]
		upper=upper[pivot]
		imod=imod[pivot]
		means=means[pivot]
		sigma=sigma[pivot, pivot]
	}
#	print(sigma)
#	print(sroot)

	e <- sapply(1:k, function(j)
	{
		x <- which(sroot[j, 1:j] != 0)
		if(length(x) != 0) { ret <- max(x) } else {warning("zero row") ; ret <- -1 }
		ret
	})
	list.e <- lapply(1:k, function(j){ which(e == j)})
#	print(list.e)
	for(i in 1:k)
	{
		if(lower[i] >= upper[i]) stop("lower bound <= upper bound")
		if(imod[i])
		{
			if(lower[i] < 0) lower[i] <- 0
			if(upper[i] < 0) stop("Imod=T with negative upper bound !")
		}
	}
#	print(rbind(lower, upper))
	
	#if(!is.loaded("tmvnlib.so")) dyn.load("tmvnlib.so")
	xx <- matrix(0, nsamp, k)
	wts <- rep(0, nsamp)
	ret <- .Fortran("rtmvnghk",
                              n      = as.integer(nsamp),
                              d      = as.integer(k),
                              means  = as.double(means),
                              sroot  = as.double(sroot),
                              a      = as.double(lower), 
                              b      = as.double(upper),
                              imod   = as.integer(imod),
                              elen   = as.integer(sapply(list.e, length)),
                              epos   = as.integer(unlist(list.e)),
                              X      = as.double(xx),
                              W = as.double(wts),
                              NAOK=TRUE, PACKAGE="tmvnsim")
	#if(is.loaded("tmvnlib.so")) dyn.unload("tmvnlib.so")
	xx <- matrix(ret$X, ncol=k, byrow=TRUE)
	wts <- ret$W
	if(any(is.na(xx)) || any(is.nan(xx)) || any(is.na(wts)) || any(is.nan(wts))) stop("NA/NaN in wts !!")  
	return(list(samp=xx[, ord], wts=wts))	
}

