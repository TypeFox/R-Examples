###rv5 with Bootstrap Censored Weibull Mixture
bootstrapCenWbMix <- function(dat, qInt=0.05, canSet=c(0.5, 0.7, 1), B=1000, iniVec=NULL,
							randSeed=NULL, conCr=1e-6, nIter=1000)
{
	##canSet is proportion
	precheck <- c(dat, canSet, qInt)
	if (any(is.na(precheck))||any(precheck<0))
		stop("Wrong Input. NO NA/Negatives allowed in the input")
	dat <- sort(dat)
	N <- length(dat)
	K <- length(canSet)
	indSet <- round(N*canSet)
	CxSet <- dat[indSet]
	
	if (is.null(iniVec))
	{	
		ttini <- genWbMixIni()
		iniVec <- rep(ttini, times=K)
	}
	else
	{
		if (length(iniVec)!=(6*K))
			stop("Wrong Input of the initial values")
	}
	
	resMx <- matrix(NA, K, 4)
	resMx[,1] <- canSet
	colnames(resMx) <- c("Censor.Prop", "Quantile.Est", "Bootstrap.S.E.", "Boostrap.RMSE")
	parmMx <- matrix(NA, K, 6)
	colnames(parmMx) <- c("Prop1", "Prop2", "Shape1", "Shape2", "Scale1", "Scale[2]")
	for (k in 1:K)
	{
		tini <- iniVec[c(((k-1)*6+1):(k*6))]
		tmix <- .Call("R2C_CWbMix", dat[c(1:indSet[k])], CxSet[k], N, indSet[k], tini, conCr, nIter)
		if (tmix[1]==0)
		{
			parmMx[k, ] <- tmix[-c(1,2)]
			tini <- tmix[-c(1,2)]
			iniVec[c(((k-1)*6+1):(k*6))] <- tini
			resMx[k, 2] <- quanWbMix.int(qInt, tmix[-c(1,2)])
		}
	}
	
	##The bootstrap part
	if (!is.null(randSeed))
		set.seed(randSeed)
	bQuanMx <- matrix(NA, B, K)
	bQuanVec<- .Call("R2C_bstpWbMix", dat, qInt, indSet, iniVec, B, conCr, nIter)
	
	
	bQuanVec[bQuanVec<0] <- NA
	bQuanMx <- matrix(bQuanVec, ncol=K, byrow=TRUE)
	
	##Processing the results
	resMx[,3] <- apply(bQuanMx, 2, sd, na.rm=T)
	empQ <- quantile(dat, qInt, type=9)
	resMx[,4] <- apply(bQuanMx, 2, mse, tq=empQ)
	
	return(list(results=resMx, parameters=parmMx, bQEst=bQuanMx))
}