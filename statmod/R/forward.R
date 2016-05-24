forward <- function(y,x,xkept=NULL,intercept=TRUE,nvar=ncol(x))
#	Forward selection for linear regression
#	30 Jan 2013
{
#	Check y
	y <- as.numeric(y)
	n <- length(y)

#	Check x
	x <- as.matrix(x)
	if(nrow(x) != n) stop("nrow of x must match length of y")

#	Check xkept
	if(!is.null(xkept)) {
		xkept <- as.matrix(xkept)
		if(nrow(xkept) != n) stop("nrow of xkept must match length of y")
	}

#	Add intercept
	if(intercept) xkept <- cbind(rep.int(1,n),xkept)
	
#	Sweep out xkept columns
	if(is.null(xkept)) {
		rank.xkept <- 0
	} else {
		QR <- qr(xkept)
		y <- qr.resid(QR,y)
		x <- qr.resid(QR,x)
		rank.xkept <- QR$rank
	}

#	Check nvar
	nvar <- min(nvar,ncol(x),n-rank.xkept)
	if(nvar <= 0) return(numeric(0))

	orderin <- rep.int(0,nvar)
	candidates <- 1:ncol(x)
	for (nin in 1:nvar) {
		if(ncol(x)==1) {
			orderin[nin] <- candidates
			break
		}

#		Standardize
		y <- y/sqrt(sum(y^2))
		x <- t(t(x)/sqrt(colSums(x^2)))

#		Next to add
		b.y.x <- crossprod(x,y)
		bestj <- which.max(abs(b.y.x))
		bestx <- x[,bestj]

#		Record and remove best covariate
		orderin[nin] <- candidates[bestj]
		candidates <- candidates[-bestj]
		x <- x[,-bestj,drop=FALSE]

#		Orthogonalize remaining wrt best covariate
		y <- y - b.y.x[bestj]*bestx
		b.x.x <- crossprod(x,bestx)
		x <- x - matrix(bestx,ncol=1) %*% matrix(b.x.x,nrow=1)
	}
	orderin
}


