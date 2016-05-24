######Desgin posterior of Eta and Beta.
######Updated on 0406, fixed the coding mistake in EtaGBeta
######Updated on Sep04, EtaGBeta part is removed to the gbase3 and merged with the NLLH function.
######Updated on Sep06, fixed a coding mistake in the B matrix part.
#######Only the function to deal with the full sequence is included.
#######name also changed, _ is used for the conditional sign.

#source("gBase3.r")

postMar.BB.Eta<- function(piMx, XMx, gpsList, rFullList=FALSE, printK=FALSE)
{
	#If piMx has one row. It is the empirical Bayesian. Only the parameter is inputed
	#Changed to version where we update (beta, etaG|s2D, s2H)
	#piMx is fixed with three column, s2H, s2D, wVec
	Q <- nrow(piMx)#Number of parameter grid.
	if (is.null(Q))
	{
		piMx <- matrix(c(piMx, 1), nrow=1)
		Q <- 1
	}
	if (is.null(gpsList$dMx))
		NoBetas <- 0
	else
		NoBetas <- ncol(gpsList$dMx) #Number of betas.
	SigSize <- NoBetas +2
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
	wVec <- piMx[, 3]
	
	##Mean preparation
	adjVec <- (Atime-Atime[1])/(Atime[T]-Atime[1])
	fVec <- Y[1] + adjVec*(Y[K] - Y[1])
	gVecGPS.cpM <- matrix(0, K, Q) ##Originally the xMeanGPS, store the Z_G beta
	if (NoBetas>0)
	{
		dMx <- as.matrix(XMx[,-c(1, 2)])
		dMxG <- gpsList$dMx
	}
	if(printK)	print("Mean Level Preparation Done")
	
	###GPS points, need to store every mean/var/cov conditioning on each phi
	etaGMean.cpM <- matrix(NA, K, Q) #Storage of the posterior mean|phi, condition on phi.
	etaGVar.cpM <- matrix(NA, K, Q) #Store the posterior var|phi
	etaGCov1.cpM <- matrix(NA, K, Q)
	betaMean.cpM <- matrix(NA, NoBetas, Q)
	if (NoBetas>0)
	{
		betaCov.cpA <- array(NA, dim=c(NoBetas, NoBetas, Q))
		betaEtaCov.cpA <- array(0, dim=c(NoBetas, K, Q))
	}
	for (i in 1:Q)
	{
		tgpost <- nllh.BB.Phi_XY(piMx[i, 1:2], gpsList, retPost=TRUE)
		etaGMean.cpM[,i] <- tgpost$etaG[,1]
		etaGVar.cpM[, i] <- tgpost$etaG[,2]
		etaGCov1.cpM[,i] <- tgpost$etaG[,3]
		if (NoBetas>0)
		{
			tBeta <- tgpost$betaMean
			gVecGPS.cpM[, i] <- dMxG %*% tBeta
			betaMean.cpM[,i] <- tBeta
			betaCov.cpA[, , i] <- tgpost$betaCov[1:NoBetas, 1:NoBetas]
			betaEtaCov.cpA[, 2:(K-1), i] <- tgpost$betaCov[1:NoBetas, -c(1:NoBetas)]
		}
	}
	rm(tgpost)
	etaGPost <- etaGMean.cpM%*% wVec
	etaGPostMeanDiffMx <- etaGMean.cpM - as.vector(etaGPost)
	etaGPostVar <- (etaGVar.cpM + etaGPostMeanDiffMx^2)%*%wVec
	etaPost <- rep(NA, T)
	etaPost[Gind] <- etaGPost 
	etaPostVar <- rep(NA, T)
	etaPostVar[Gind] <- etaGPostVar
	if(printK)	print("GPS points Done")
	if (rFullList)
	{
		etaPostF <- matrix(NA, T, Q)
		etaPostVarF <- matrix(NA, T, Q)
		etaPostF[Gind, ] <-  etaGMean.cpM
		etaPostVarF[Gind, ] <- etaGVar.cpM
	}
	############################################################################
	###Calculation begin for the non-GPS points.
	for (k in 1:(K-1))
	{
		if (printK)	print(paste("Period ", k))
		tSind <- Gind[k]
		tEind <- Gind[(k+1)]
		if (tEind>(tSind+1))
		{
			#Common for all the Q pi choices.
			tInd <- (tSind+1):(tEind-1)
			tGind <- c(tSind, tEind)
			tT <- length(tInd)
			tTime <-Atime[tInd]	
			tAdjVec <- (tTime - Gtime[k])/(Gtime[k+1] - Gtime[k])
			tEtaOrgCov <- covBBridge.diag(c(Gtime[k], tTime, Gtime[k+1]))
						#The original covariance matrix of the Brownian Bridge.
			
			tEtaMean.cGcX <- matrix(NA, tT, Q) #Each column conditioning on phi
			tEtaVar.cGcX <- matrix(NA, tT, Q)
			
			for (i in 1:Q)
			{
				tS2H <- piMx[i, 1]
				tS2D <- piMx[i, 2]
				
				tAmx <- cbind(1-tAdjVec, tAdjVec)
				tFvec <- tAmx %*% c(etaGMean.cpM[k, i], etaGMean.cpM[k+1, i])
					#Following the notation in the proposal
				tBeta <- betaMean.cpM[,i]
				tgtid <- rep(0, tT)
				if (NoBetas >0)
					tgtid <- matrix(dMx[tInd, ], ncol=NoBetas) %*% tBeta -
								tAmx %*% gVecGPS.cpM[c(k, k+1), i]
				tXtid <- tAmx %*% X[tGind]
				
				tSnRatio <- tS2H/(tS2D+tS2H) #signal noise ratio
				tEtaMean.cGcX[,i] <- tFvec + tSnRatio*(X[tInd] - tgtid- tXtid)
				
				#Variance From the conditional distribution
				tEtaDiagCov.cGcX <- tEtaOrgCov*tS2H*(1-tSnRatio)
				#Variance From marginalize the etaG, beta
				tBmx <- tAmx
				tSig <- matrix(0, SigSize, SigSize)
				tSig[1:2, 1:2] <- matrix(c(etaGVar.cpM[k, i], etaGCov1.cpM[k, i], 
								etaGCov1.cpM[k, i],etaGVar.cpM[(k+1), i]), ncol=2)
				tZ <- NULL
				if(NoBetas >0)
				{
					tSig[3:SigSize, 3:SigSize] <- betaCov.cpA[, , i]
					tSig[3:SigSize, 1:2] <- betaEtaCov.cpA[,c(k, k+1), i]
					tSig[1:2, 3:SigSize] <- t(tSig[3:SigSize, 1:2])
					
					tZ <- tSnRatio*matrix(dMx[tInd, ]- tAmx %*% dMxG[c(k, k+1), ], ncol=NoBetas)
					tBmx <- cbind(tBmx, tZ)
				}
				if (k==1||k==(K-1))
				{ ###Forgot rho part here? Sep4.....  
					tSigEig <- eigen(tSig)
					#print(tSigEig)
					tval <- tSigEig$val
					tval[tval<0] <-0
					tSigCM <- tSigEig$vectors %*% diag(sqrt(tval))
					tBV <- tBmx %*% tSigCM
				}
				else
				{
					tSigCM <- chol(tSig)
					tBV <- tBmx %*% t(tSigCM)
				}
				tEtaVar.cGcX[, i] <- tEtaDiagCov.cGcX + rowSums(tBV^2)
				rm(tBV)
			}
			if (rFullList)
			{
				etaPostF[tInd, ] <-  tEtaMean.cGcX
				etaPostVarF[tInd, ] <- tEtaVar.cGcX
			}
			tEtaPost<- tEtaMean.cGcX %*% wVec
			etaPost[tInd]  <- tEtaPost
			tEtaMeanDiffMx <- tEtaMean.cGcX - as.vector(tEtaPost)
			etaPostVar[tInd] <- (tEtaVar.cGcX + tEtaMeanDiffMx^2)%*%wVec
						
			rm(tSind, tEind, tInd, tGind, tT, tTime, tAdjVec, tEtaOrgCov,
			tEtaMean.cGcX, tEtaVar.cGcX, tS2H, tS2D, tAmx, tFvec, tBeta, tgtid,
			tXtid,tSnRatio,tBmx, tSig, tZ, tEtaPost, tEtaMeanDiffMx)
		}
	}
	if (rFullList)
	{
		return(list(Mean.cpM=etaPostF, Var.cpm=etaPostVarF, 
				wVec=wVec, Margin=cbind(Mean=etaPost, Var=etaPostVar)))
	}
	else
		return(cbind(Mean=etaPost, Var=etaPostVar))
}
    