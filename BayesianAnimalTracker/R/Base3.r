####Base Function. The covariance matrix, random number generations..
####Time starts at 1.
####uses both for the conditional approach and direct approach.
#library(Matrix)
####Modified on Sep 3rd.


offDiag1 <- function(m)
{
	k <- nrow(m)
	if(k>2)
		return(diag(m[,-1]))
	if(k==1)
		return(0)
	if(k==2)
		return(m[1,2])
}

rBBridge <- function(Tind, A=0, B=0, s2H=1)
{
	#Function to simulate a Brown bridge the starting point X[1]=A!!!!!
	#X(T)=B, the base variance is s2H, variance of the incremental part in the Brownian Motion
	#Tind is the time index. begining at 1.
	K <- length(Tind)
	tstep <- Tind[-1] - Tind[-K]
	sdvec <- sqrt(s2H*tstep)
	bm <- cumsum(c(0, rnorm(K-1, sd=sdvec)))
	adjVec <- (Tind - Tind[1])/(Tind[K] - Tind[1])
	bb <- bm - bm[K]* adjVec + A + (B-A)*adjVec
	return(bb)
}

rBMotion <- function(Tind, s2D=1)
{
	#Function to simulate a Brown motion, starting at zero
	#Tind is the time index. begining at 1.
	#innovation variance s2H*tLag
	K <- length(Tind)
	tstep <- Tind[-1] - Tind[-K]
	sdvec <- sqrt(s2D*tstep)
	bm <- cumsum(c(0, rnorm(K-1, sd=sdvec)))
	return(bm)
}

covBBridge <- function(Tind)
{
	#The covariance matrix of a Brownian Bridge.
	#Return a matrix of T-2, time starts at 1.
	K <- length(Tind)
	m <- matrix(0, K, K)
	totT <- Tind[K] - Tind[1]
	for (i in 2:(K-1))
	{
		for (j in i:(K-1))
		{
			m[i, j] <- (Tind[i]-Tind[1])*(Tind[K]-Tind[j])/totT
			m[j, i] <- m[i, j]
		}
	}
	return(m[2:(K-1), 2:(K-1)])
}

corrAR1 <- function(Tind, rho)
{
	#A little difference from covBBridge, returns a matrix of K-1
	#input Tind is the same. The first point T=1 is discarded.
	K <- length(Tind)
	m <- diag(K)
	for (i in 1:(K-1))
	{
		for (j in (i+1):K)
		{
			m[i, j] <- m[j, i] <- rho^(abs(Tind[i]-Tind[j]))
		}
	}
	return(m[2:K, 2:K])
}

covBMotion <- function(Aind)
{
	T <- length(Aind)
	m <- matrix(NA, T, T)
	for(i in 1:T)
	{
		m[1:i,i] <- Aind[1:i]-Aind[1]
		m[-c(1:i), i] <- Aind[i]-Aind[1]
	}
	return(m[2:T, 2:T])
}

covBBridge.diag <- function(Tind)
{
	#The diagonal element of the covariance matrix of a Brownian Bridge.
	#Return a vector or length T
	K <- length(Tind)
	totT <- Tind[K] - Tind[1]
	m <- (Tind-Tind[1])*(Tind[K] - Tind)/totT
	return(m[-c(1,K)])
}

###The ad-hoc bias correction method from the biology literature.
adhocEta <- function(XMx, gpsList)
{
	T <- nrow(XMx)	
	K <- gpsList$K
	Y <- gpsList$Y
	X <- XMx[,1]
	Gtime <- gpsList$Gtime
	Atime <- XMx[,2]
	if ((Atime[1]!=Gtime[1])||(Atime[T]!=Gtime[K]))
		stop("Input Time Error")
	Aind <-1:T
	#Gind <- Aind[(Atime%in%Gtime)]
	Gind <- gpsList$Gind
	##Mean preparation
	etaBC <- rep(NA, T)
	etaBC[Gind] <- Y
	for (k in 1:(K-1))
	{
		tSind <- Gind[k]
		tEind <- Gind[(k+1)]
		if (tEind>(tSind+1))
		{
			tInd <- (tSind+1):(tEind-1)
			tBiasCoc <- X[tEind]- X[tSind]
			tT <- length(tInd)
			tTime <-Atime[tInd]	
			tAdjVec <- (tTime - Gtime[k])/(Gtime[k+1] - Gtime[k])
			tfVec <- Y[k]*(1-tAdjVec) + Y[(k+1)]*tAdjVec
			tXvec <- (X[tInd] - X[tSind]) - tAdjVec*tBiasCoc
			etaBC[tInd] <-tfVec + tXvec
		}
	}
	return(etaBC)
}


xvAnn <- function(mat)
{
	if (ncol(mat)<3)
		return(mean((mat[,1]-mat[,2])^2))
	else
	{
		cvrmse <- mean((mat[,1]-mat[,2])^2)
		cat("XV-rmse \n")
		print(cvrmse)
		ubound <- mat[,2] + sqrt(mat[,3])*1.96
		lbound <- mat[,2] - sqrt(mat[,3])*1.96
		coverageNum <- sum((mat[,1] <=ubound)*(mat[,1] >=lbound)) 
		cat("Coverage Percentage \n")
		per <- coverageNum/nrow(mat)
		print(per)
		return(list(cvrmse=cvrmse, coverageNum=coverageNum, per=per, mat=cbind(mat, lbound, lbound) ))
	}	
}


linInter <- function(Ntime, Stime, Etime, Sloc, Eloc)
{
	ydif <- Eloc-Sloc
	xdif <- Etime-Stime
	return(Sloc + (Ntime- Stime)*ydif/xdif)
}


xvAnn2 <- function(mat)
{
	##CV_RMSE
	if (ncol(mat)<3)
		return(sqrt(mean((mat[,1]-mat[,2])^2)))
	else
	{
		cvrmse <- sqrt(mean((mat[,1]-mat[,2])^2))
		cat("XV-rmse \n")
		print(cvrmse)
		ubound <- mat[,2] + sqrt(mat[,3])*1.96
		lbound <- mat[,2] - sqrt(mat[,3])*1.96
		coverageNum <- sum((mat[,1] <=ubound)*(mat[,1] >=lbound)) 
		cat("Coverage Percentage \n")
		per <- coverageNum/nrow(mat)
		print(per)
		return(list(cvrmse=cvrmse, coverageNum=coverageNum, per=per, mat=cbind(mat, lbound, lbound) ))
	}	
}