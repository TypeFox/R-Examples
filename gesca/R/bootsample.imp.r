bootsample.imp <- function (z0, case_index, nvar, nobs_g, ng, b, nobs_tot, moption, missingvalue)
{
	#---------------------------------------
	# output:
 	# 	Z = block diagonal matrix of "normalized" data after mean substitution
 	# 	z0_meanimp = nobs_tot by nvar matrix of data after mean substitution
 	# 	rawz0 = block diagonal matrix of original data with missings
 	# Last revised July 28, 2015
 	# used only for original data (b == 0)
	#---------------------------------------

	if (moption == 3) rawz0 <- matrix(0, nobs_tot, ng*nvar)
	
	Z <- matrix(0, nobs_tot, ng*nvar)
	z0_meanimp <- matrix(0, nobs_tot, nvar)
	kk <- 0
	for (g in 1:ng) {
	    k <- kk + 1
		kk <- kk + nvar
		if ( b == 0 ) {
			bz01 <- z0[case_index[g,1]:case_index[g,2],]
		} else {
	       rb <- runif(nobs_g[g])*nobs_g[g]
		   rrb <- ceiling(rb)
		   bz01 <- matrix(0, nobs_g[g], nvar)
		   z0_g <- z0[case_index[g,1]:case_index[g,2],]
		   for (i in 1:nobs_g[g]) { bz01[i,] <- z0_g[rrb[i],] }
		}
		
		if ( moption == 3 ) {
			rawz0[case_index[g,1]:case_index[g,2],k:kk] <- bz01
		} else {
			rawz0 <- c()
		}
	
		# Columnwise mean substitution for each group
		for (j in 1:nvar) {
			bz01_col <- bz01[,j]
			indcol <- which(bz01_col == missingvalue)
			nnmissingobs <- nobs_g[g] - length(indcol)
			bz01_col[indcol] <- 0
			mean_col <- sum(bz01_col)*(1/nnmissingobs)
			bz01_col[indcol] <- mean_col
			bz01[,j] <- bz01_col
		}
		z0_meanimp[case_index[g,1]:case_index[g,2],] <- bz01 
		ctz <- bz01 - ( matrix(1,nobs_g[g],1) %*% apply(bz01,2,mean) )
		covz <- (t(ctz)%*%ctz)/nobs_g[g]
		Dstz <- sqrt(diag(diag(covz)))
		bz0 <- t(solve(t(Dstz),t(ctz))) #  bz0 = ctz/Dstz in matlab
		bz0 <- bz0/sqrt(nobs_g[g]) # data normalization after mean substitution
		Z[case_index[g,1]:case_index[g,2],k:kk] <- bz0
	
		rm(bz01)
		rm(bz0)
	}
	
	output.bootsample.imp <- list(Z = Z, z0_meanimp = z0_meanimp, rawz0 = rawz0)
	output.bootsample.imp
}