####Censored Weibull MLE with bootstrap
mse <- function(vec, tq) return(sqrt(mean((vec-tq)^2, na.rm=T)))

bootstrapCMLE <- function(dat, qInt=0.05, canSet=seq(0.1,0.5, by=0.1), B=5000, 
							randSeed=NULL, conCr=1e-9, nIter=1000)
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
	
	resMx <- matrix(NA, K, 6)
	resMx[,1] <- canSet
	colnames(resMx) <- c("Censor.Prop", "Shape.Est", "Scale.Est", "Quantile.Est", "Bootstrap.SE", "Boostrap.RMSE")
	for (k in 1:K)
	{
		tcmle <- .Call("R2C_cenWeibullMLE", dat[c(1:indSet[k])], CxSet[k], N, indSet[k], conCr, nIter)
		if (tcmle[1]==0)
		{
			resMx[k, 2:3] <- tcmle[-1]
			resMx[k, 4] <- qweibull(qInt, tcmle[2], tcmle[3])
		}
	}
	##The bootstrap part
	if (!is.null(randSeed))
		set.seed(randSeed)
	bQuanMx <- matrix(NA, B, K)
	bQuanVec<-.Call("R2C_boostrapCMLE", dat, qInt, indSet, B, conCr, nIter)
	bQuanVec[bQuanVec<0] <- NA
	bQuanMx <- matrix(bQuanVec, ncol=K, byrow=TRUE)
	
	##Processing the results
	resMx[,5] <- apply(bQuanMx, 2, sd, na.rm=T)
	empQ <- quantile(dat, qInt, type=9)
	resMx[,6] <- apply(bQuanMx, 2, mse, tq=empQ)
	
	return(list(results=resMx, bQEst=bQuanMx))
}
