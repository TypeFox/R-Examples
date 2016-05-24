evalOC <-
function(betaTruth,
									 betaHat,
									 seHat,
									 waldTest,
									 keep,
									 alpha=0.05,
									 referent=2,
									 resNames=NULL)
{	
	## Useful objects
	##
	nDesigns <- dim(betaHat)[2]
	p        <- dim(betaHat)[3]
	betaZero <- betaTruth
	betaZero[which(betaZero == 0)] <- NA
	
	##
	betaMean       <- matrix(NA, nrow=nDesigns, ncol=p, dimnames=resNames)
	betaMeanBias   <- matrix(NA, nrow=nDesigns, ncol=p, dimnames=resNames)
	betaMeanPB     <- matrix(NA, nrow=nDesigns, ncol=p, dimnames=resNames)
	betaMedian     <- matrix(NA, nrow=nDesigns, ncol=p, dimnames=resNames)
	betaMedianBias <- matrix(NA, nrow=nDesigns, ncol=p, dimnames=resNames)
	betaMedianPB   <- matrix(NA, nrow=nDesigns, ncol=p, dimnames=resNames)
	betaSD         <- matrix(NA, nrow=nDesigns, ncol=p, dimnames=resNames)
	betaMSE        <- matrix(NA, nrow=nDesigns, ncol=p, dimnames=resNames)
	seMean         <- matrix(NA, nrow=nDesigns, ncol=p, dimnames=resNames)
	seRatio        <- matrix(NA, nrow=nDesigns, ncol=p, dimnames=resNames)
	betaCP         <- matrix(NA, nrow=nDesigns, ncol=p, dimnames=resNames)
	betaPower      <- matrix(NA, nrow=nDesigns, ncol=p, dimnames=resNames)
	##
	for(i in 1:nDesigns)
	{
		##
		if(sum(keep[,i]) > 0)
		{
			##
			tempB  <- betaHat[keep[,i],i,]
			tempSE <- seHat[keep[,i],i,]
			tempP  <- waldTest[keep[,i],i,]
			##
			betaMat <- matrix(betaTruth, nrow=nrow(tempB), ncol=p, byrow=TRUE)
		
  		##
  		betaMean[i,]       <- apply(tempB, 2, mean)
  		betaMeanBias[i,]   <- betaMean[i,] - betaTruth
			betaMeanPB[i,]     <- (betaMean[i,] - betaZero) / betaZero * 100
  		betaMedian[i,]     <- apply(tempB, 2, median)
  		betaMedianBias[i,] <- betaMedian[i,] - betaTruth
			betaMedianPB[i,]   <- (betaMedian[i,] - betaZero) / betaZero * 100
			betaSD[i,]         <- apply(tempB, 2, sd)
			betaMSE[i,]        <- betaMeanBias[i,]^2 + apply(tempB, 2, var)
			seMean[i,]         <- apply(tempSE, 2, mean)
			seRatio[i,]        <- seMean[i,] / betaSD[i,] * 100
			ciL                <- tempB + (qnorm(alpha/2) * tempSE)
			ciU                <- tempB - (qnorm(alpha/2) * tempSE)
			betaCP[i,]         <- apply(((ciL < betaMat) & (betaMat < ciU)), 2, mean) * 100
			betaPower[i,]      <- apply(tempP, 2, mean) * 100
		}
	}

  ## Relative Uncertainty
  betaRU <- betaSD / matrix(betaSD[referent,], nrow=nDesigns, ncol=p, byrow=TRUE) * 100
	dimnames(betaRU) <- resNames
	
	##
	value <- list(betaMean=betaMean,
								betaMeanBias=betaMeanBias,
								betaMeanPB=betaMeanPB,
								betaMedian=betaMedian,
								betaMedianBias=betaMedianBias,
								betaMedianPB=betaMedianPB,
								betaSD=betaSD,
								betaMSE=betaMSE,
								seMean=seMean,
								seRatio=seRatio,
								betaCP=betaCP,
								betaPower=betaPower,
								betaRU=betaRU)
  ##
  return(value)
}
