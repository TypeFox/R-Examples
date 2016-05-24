####Modified on Nov11-2014. 
####Simplified for package.
dataSim <- function(T, K, s2H, s2D, s2G=0.01, 
		gind=NULL, betaVec=NULL, dMx=NULL, A=0, B=0, scale=TRUE)
{
	Aind <- 1:T
	Atime <- Aind
	Q <- length(betaVec)
	#dMx can be inputed as a T by Q matrix.
	
	if (Q>0)
	{
		if (is.null(dMx))
		{
			dMx <- matrix(NA, nrow=T, ncol=Q)
			for (j in 1:Q)
			{
				tx <- Atime^(j-1)
				if (j >1&scale)
					tx <- scale(tx)
				dMx[, j] <- tx
			}
		}
		gVec <- dMx %*% betaVec
	}
	else
	{
		dMx <- NULL
		gVec <- rep(0, T)
	}
	
	if (is.null(gind))
		Gind <- c(1, sort(sample(2:(T-1), K-2, replace=FALSE)), T)
	else
		Gind <- c(1, gind, T)
	Gtime <- Gind
	
	eta <- rBBridge(Atime, A, B, s2H=s2H)
	y <- eta[Gind]
	y[2:(K-1)] <- y[2:(K-1)] + rnorm(K-2, sd=sqrt(s2G))
	
	xi <- rBMotion(Atime, s2D) + gVec
	x <- xi + eta
	return(list(eta=eta, Y=y, Ytime=Gtime, X=x)) 
}


as.dataList <- function(X, Y, Ytime, Xtime=NULL, s2G, timeUnit=1, dUnit=1,
		 dMx=NULL, betaOrder=1, scale=TRUE)
{
	X <- X/dUnit
	Y <- Y/dUnit
	T <- length(X)
	K <- length(Y)
	
	if (is.null(Xtime))
		Xtime <- 1:T
	Xind <- 1:T	
	Yind <- Xind[Xtime %in% Ytime]	
	
	if (any(class(Xtime) %in% c("POSIXlt", "POSIXt")))
	{	
		Xtime <- as.numeric(Xtime - Xtime[1])+1
		Ytime <- Xtime[Yind]
	}
	else
	{
		Xtime <- Xind
		Ytime <- Yind
	}
	Xtime <- Xtime/timeUnit
	Ytime <- Ytime/timeUnit
	
	R0 <- covBBridge(Ytime)
	adjVec <- (Ytime - Ytime[1])/(Ytime[K] - Ytime[1])
	fVec <- Y[1] + adjVec*(Y[K] - Y[1])
	
	pMxBB <- solve(R0)
	pMxBM <- solve(covBMotion(Ytime))
	
	if (is.null(dMx))
	{
		if (betaOrder==0)
			dMx <- NULL
		else
		{
			dMx <- matrix(NA, nrow=T, ncol=betaOrder)
			for (j in 1:betaOrder)
			{	
				tx <- Xtime^(j-1)
				if (j >1&scale)
					tx <- scale(tx)
				dMx[, j] <- tx
			}
		}
	}
	else
	{
		dMx <- matrix(dMx, nrow=T)
		betaOrder=ncol(dMx)
	}
	glist <- list(X=X[Yind], Y=Y, Gtime=Ytime, Gind=Yind, s2G=s2G, K=K, 
				R0=R0, fVec=fVec, dMx=NULL, pMxBB=pMxBB, pMxBM=pMxBM)
	if (betaOrder>0)
		glist$dMx <-as.matrix(dMx[Yind, ], nrow=K, ncol=betaOrder)
	
	XMx <- cbind(X, Xtime, dMx)
	return(list(XMx=XMx, glist=glist))
}