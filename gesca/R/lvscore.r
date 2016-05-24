lvscore <- function (x, W, option = 1)
{
	#---------------------------------------------------
	# calculate latent scores
	# two approaches: option = 1 (Fornell); = 2 (Tenenhaus)
	# Last revised Aug 27, 2015
	#---------------------------------------------------

	nobs <- nrow(x)
	nlv <- ncol(W)
	ctz <- x - matrix(1,nobs,1)%*%apply(x,2,mean)
	covz <- (t(ctz)%*%ctz)/nobs
	Dstz <- sqrt(diag(diag(covz)))
	UstdW <- solve(Dstz,W)
	if ( option == 1 ) {
		sumUstdW <- apply(UstdW,2,sum)
		rUstdW <- UstdW 			# relative unstandardized weights
		for (j in 1:nlv) {
			rUstdW[,j] <- UstdW[,j]/sumUstdW[j]
		}
		lvscore <- x%*%rUstdW 		# Fornell's approach
	} else if ( option == 2 ) {
		lvscore <- x%*%UstdW 		# Tenenhaus's approach
	}
	
	lvscore
}