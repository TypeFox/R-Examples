modelfit.mg <- function (Z0, W, T, nvar, nlv, ng, case_index)
{
	# MODEL FIT MEASURES
	gfi_1 <- 0
	gfi_2 <- 0
	srmr_1 <- 0
	kk <- 0
	ss <- 0
	COR_RES <- c() # correlation residual matrix
	for (g in 1:ng) {
		k <- kk + 1
		kk <- kk + nvar
		s <- ss + 1
		ss  <- ss + nlv
		
		zz <- Z0[case_index[g,1]:case_index[g,2],k:kk]
		w <- W[k:kk, s:ss]
		v <- cbind(diag(nvar), w)
		t <- T[s:ss,]
		omega <- v - w%*%t
		ee <- zz%*%omega
		
		samcov <- cov(zz)	# sample covariances for each group
		samcorr <- cor(zz)
		tp_precov <- solve(omega%*%t(omega),omega%*%diag(apply(ee,2,var))%*%t(omega))
		precov <- t(solve(t(omega%*%t(omega)),t(tp_precov)))
		COV_RES <- samcov - precov
		prerij <- precov
		for (i in 1:nvar) {
			for (j in 1:nvar) { prerij[i,j] <- precov[i,j]/sqrt(precov[i,i]*precov[j,j]) }
		}
		srmr <- 0
		for (i in 1:nvar) {
			for (j in 1:nvar) {
				if (j > i) {
					corr_residual <- (samcorr[i,j] - prerij[i,j])^2
					srmr <- srmr + corr_residual
				}
			}
		}
		srmr_1 <- srmr_1 + srmr
		gfi_1 <- gfi_1 + sum(diag(COV_RES^2))
		gfi_2 <- gfi_2 + sum(diag(samcov^2))
		COR_RES <- rbind(COR_RES,samcorr - prerij)
	}
	
	nvar_tot <- ng*nvar
	srmr_2 <- nvar_tot*(nvar_tot+1)/2
	SRMR <- sqrt(srmr_1/srmr_2)				# Standardized root mean square residual
	GFI <- 1 - (gfi_1/gfi_2)				# GFI-ULS
	
	output.modelfit.mg <- list(GFI = GFI, SRMR = SRMR, COR_RES = COR_RES)
	output.modelfit.mg
}