####0903-Evaluate [phi|X_G, Y] as in the proposal
#source("Base3.r")

nllh.BB.Phi_XY<- function(parm, gpsList, logS2=FALSE, retPost=FALSE)
{
	#This is the combination of the gBase2 nllh function and fBase2 post.EtaGBeta function.
	#It is implemented by marginalize all the beta, eta and then evaluate phi. 
	#The conditional posterior of beta, etaG could also be returned including their covariance lag 0,1 
	#Including the variance and covariance of lag 1 in the return.
	if (logS2)
	{
		s2H <- exp(parm[1])
		s2D <- exp(parm[2])
	}
	else
	{
		s2H <- parm[1]
		s2D <- parm[2]
		if (any(parm[1:2]<0))
			return(1E30)
	}
	Y <- gpsList$Y
	X <- gpsList$X 
	K <- gpsList$K
	
	Gtime <- gpsList$Gtime
	pMx.BB <- gpsList$pMxBB/s2H
	pMx.WN <- diag(K-2)/gpsList$s2G #White noise
	pMx.X <- gpsList$pMxBM/s2D
	Xstar <- X[-1]
	Xstar[K-1] <- Xstar[K-1] - Y[K]
	
	##Eta star  c(\beta, etaG) length Q + K-2
	if (is.null(gpsList$dMx))
		Q <-0
	else
		Q <- ncol(gpsList$dMx)
	
	G <- cbind(matrix(0, nrow=K-2, ncol=Q), diag(K-2))
	D <- cbind(gpsList$dMx[-1, ], rbind(diag(K-2), rep(0, K-2)))
	
	gpbb <- t(G)%*%pMx.BB
	gpwn <- t(G)%*%pMx.WN
	dpx <-  t(D)%*% pMx.X
	fg <-  gpsList$fVec[-c(1,K)]
	yg <-  Y[-c(1, K)]
	
	M1 <- gpbb%*%G + gpwn%*%G +dpx%*%D 
	M2 <- gpbb %*%fg + gpwn %*%yg + dpx %*% Xstar
	m31 <- fg%*%pMx.BB%*%fg + t(yg)%*%pMx.WN%*%yg + t(Xstar)%*% pMx.X%*%(Xstar)
	if(retPost)
	{
		M1Inv <- solve(M1)
		m3 <- m31 - t(M2) %*% M1Inv%*% M2
	}
	else
		m3 <- m31 - t(M2) %*% solve(M1, M2)
	
	det1 <- -determinant(pMx.BB, logarithm = TRUE)$modulus/2
	det2 <- -determinant(pMx.WN, logarithm = TRUE)$modulus/2
	det3 <- -determinant(pMx.X, logarithm = TRUE)$modulus/2
	det4 <- determinant(M1, logarithm = TRUE)$modulus/2
	nllh <- det1 + det2 + det3 +det4
	attr(nllh, "logarithm") <- NULL
	nllh <- nllh + m3/2
	if(retPost)
	{
		etaSPost <- M1Inv %*% M2
		if (Q>0)
		{
			etaGPost <- etaSPost[-c(1:Q)]
			etaGPostCov <- M1Inv[-c(1:Q), -c(1:Q)]
			etaGPostVar <- diag(etaGPostCov)
			etaGPostCov1 <- offDiag1(etaGPostCov)
			betaCov <- matrix(M1Inv[1:Q, ], nrow=Q, ncol=Q+K-2)
			betaMean <- etaSPost[1:Q]
		}
		else
		{
			etaGPost <- etaSPost
			etaGPostCov <- M1Inv
			etaGPostVar <- diag(M1Inv)
			etaGPostCov1 <- offDiag1(M1Inv)
			betaCov<- NULL
			betaMean <-NULL
		}
		etaGRes <- cbind(Mean=c(Y[1], as.vector(etaGPost), Y[K]), 
						Var=c(0, etaGPostVar, 0), 
						CovLag1=c(0, etaGPostCov1, 0, 0))

		return(list(nllh=c(nllh), etaG = etaGRes, betaMean=betaMean, betaCov=betaCov))
	}
	else
		return(c(nllh))
}

zSearch <-function(gnlm, glist, zStepSize=1, logPiTol=3, logPiSubTol=5)
{
	###Method adopted from INLA. 
	###Search in the direction of eigen value
	###gnlm, the nlm object of minmizing nllh.BB.marPhi logS2=TRUE
	###glist the GPS observation list
	###zStepSize, control how fine the search is
	###logPiTol, the difference in log deviance to stop search.
	###logPiSubTol, to which to abandon the points with too small LLH
	###S2 into S3 on Sep 4th.
	invHes <- solve(gnlm$hessian)
	bnllh <- gnlm$minimum
	sigEig <- eigen(invHes)
	sigEigV <- sigEig$vector
	sigEigD <- sqrt(sigEig$value)
	thetaStar <- gnlm$estimate
	P <- length(thetaStar) #P=2 for this version.
	zList <- vector("list", P)
	
	for (i in 1:P)
	{
		#print(i)
		#postive side 
		logDev <- 0
		tz <-0
		tpZVec <- NULL
		while(logDev <logPiTol)
		{
			tz <-tz + zStepSize
			tparm <- thetaStar+ tz*sigEigD[i]*sigEigV[, i]
			tnllh <- nllh.BB.Phi_XY(tparm, glist, logS2=TRUE)
			tpZVec <- c(tpZVec, tz)
			logDev <- tnllh-bnllh
		}
		#negative side search
		logDev <- 0
		tz <-0
		tnZVec <- NULL
		while(logDev <logPiTol)
		{
			tz <-tz - zStepSize
			tparm <- thetaStar+ tz*sigEigD[i]*sigEigV[, i]
			tnllh <- nllh.BB.Phi_XY(tparm, glist, logS2=TRUE)
			tnZVec <- c(tnZVec, tz)
			logDev <- tnllh-bnllh
		}
		zList[[i]] <- c(rev(tnZVec), 0, tpZVec)
	}
	zMx <- as.matrix(expand.grid(zList))
	Q <- nrow(zMx)
	res <- matrix(NA, nrow=Q, ncol=(P+1))
	for (i in 1:Q)
	{
		tparm <- thetaStar + sigEigV %*% (zMx[i, ]*sigEigD)
		res[i, 1:P] <- tparm
		res[i, P+1] <- nllh.BB.Phi_XY(tparm, glist, logS2=TRUE)
	}
	if (logPiSubTol>0)
	{
		res <- res[ ((res[,P+1]-bnllh)<=logPiSubTol), ]
	}
	#Standardize the [phi|X, Y]
	res[,P+1] <- exp(-res[,P+1])
	res[, P+1] <- res[, P+1]/sum(res[, P+1])
	#Put the two log parameter back
	res[,1] <- exp(res[,1])
	res[,2] <- exp(res[,2])
	return(res)
}
