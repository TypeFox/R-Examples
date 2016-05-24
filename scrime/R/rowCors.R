`rowCors` <-
function(X, y, trendStat=FALSE, use.n=NULL){
	if(!is.matrix(X) || !is.numeric(X))
		stop("X must be a numeric matrix.")
	if(!is.vector(y) || !is.numeric(y))
		stop("y must be a numeric vector.")
	if(ncol(X)!=length(y))
		stop("The length of y must be equal to the number of columns of X.")
	if(any(is.na(y)))
		stop("There are missing values in y.")
	if(trendStat){
		n.cat <- length(unique(y))
		if(n.cat<2)
			stop("y must consist of at least two levels.")
		if(n.cat>10)
			stop("y consists of more than 10 levels.")
		if(is.null(use.n))
			use.n <- n.cat==2
	}
	n <- rowSums(!is.na(X)) - 1 
	if(any(n<=0))
		stop("Each row of X must contain at least two non-missing values.")
	y <- (y - mean(y)) / sd(y)
	X <- X - rowMeans(X, na.rm=TRUE)
	rSd <- rowSums(X*X, na.rm=TRUE) / n
	X <- X / sqrt(rSd)
	#X[is.na(X)] <- 0
	r <- as.vector(X%*%y)
	r <- r/n
	if(!trendStat)
		return(r)
	if(use.n)
		n <- n + 1
	r*r*n
}

