`rowScales` <-
function(X, add.stats=FALSE){
	n <- rowSums(!is.na(X)) - 1
	if(any(n<=0))
		stop("Each row of X must contain at least two non-missing values.")
	rM <- rowMeans(X, na.rm=TRUE)
	X <- X - rM
	rSd <- sqrt(rowSums(X*X, na.rm=TRUE) / n)
	X <- X/rSd
	if(!add.stats)
		return(X)
	structure(list(X=X, means=rM, sds= rSd))
}

