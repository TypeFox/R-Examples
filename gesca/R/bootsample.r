bootsample <- function (z0, case_index, nvar, nobs_g, ng, b, nobs_tot)
{
	#---------------------------------------
	# only case-level data are considered
	# all bootstrap samples are normalized (MS = 0 in old codes)
	# Heungsun Hwang, Sunmee Kim
	# Last revised Aug 24, 2015
	#---------------------------------------

	Z <- matrix(0, nobs_tot, ng*nvar)
	kk <- 0
	for (g in 1:ng) {
		k <- kk + 1
		kk <- kk + nvar
		if ( b == 0) {
			bz01 <- z0[case_index[g,1]:case_index[g,2],]
		} else {
			rb <- runif(nobs_g[g])*nobs_g[g]
			rrb <- ceiling(rb)
			bz01 <- matrix(0, nobs_g[g], nvar)
			z0_g <- z0[case_index[g,1]:case_index[g,2],]
			for (i in 1:nobs_g[g]) { bz01[i,] <- z0_g[rrb[i],] }
		}
		
		# Groupwise normalization
		ctz <- bz01 - ( matrix(1,nobs_g[g],1) %*% apply(bz01,2,mean) )
		covz <- (t(ctz)%*%ctz)/nobs_g[g]
		Dstz <- sqrt(diag(diag(covz)))
		bz0 <- t(solve(t(Dstz),t(ctz)))
		bz0 <- bz0/sqrt(nobs_g[g]) # data normalization
		Z[case_index[g,1]:case_index[g,2],k:kk] <- bz0
		
	rm(bz01)
	rm(bz0)
	}
		
	output.bootsample <- list(Z = Z)
	output.bootsample
}