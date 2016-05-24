gsca.mg.ho2 <- function (z0, group_var, W001, W002, C00, B00,
				loadtype1 = matrix(1,1,ncol(W001)), loadtype2 = matrix(1,1,ncol(W002)),
				nbt = 100, itmax = 100, ceps = 0.00001, moption = 0, missingvalue = NULL)
{
	#---------------------------------------
	# Last revised September 3, 2015
	# GSCA with second-order latents (multiple-group analysis)
	# moption = 0 (no missing), 1 (listwise deletion), 2 (mean substitution), 3 (least squares imputation)
	#---------------------------------------
	
	# listwise deletion
	if ( moption == 1 ) {
		nrow <- nrow(z0)
		row_index <- rep(1,nrow)		
		for (i in 1:nrow) {
			if ( sum(which(z0[i,] == missingvalue)) > 0 ) { row_index[i] = 0 }
		}
		rindex <- which(row_index != 0)
		z0 <- z0[rindex,]
		group_var <- group_var[rindex]
	} 

	# number of groups and cases per group
	nobs_tot <- nrow(z0)
	nvar <- ncol(z0)	
	ng <- length(unique(group_var))
	nobs_g <- matrix(,1,ng)
	for (j in 1:ng) {
		nobs_g[,j] <- as.numeric(table(group_var))[j]
	}

	# case numbers for each group
	case_index <- matrix(,ng,2)
	kk <- 0
	for (j in 1:ng) {
		k <- kk + 1
		kk <- kk + nobs_g[j]
		case_index[j,1] = k
		case_index[j,2] = kk
	}
	
	# ---------------------------------------
	# model specification for a single group
	# ---------------------------------------
	
	nlv1 <- length(loadtype1)         # num of 1st-order latents
	nlv2 <- length(loadtype2)         # num of 2nd-order latents
	nlv <- nlv1 + nlv2                # total number of latents
	ntv <- nvar + nlv                 # number of indicators and latents
	num_par_W002 <- length(W002[!W002 == 0])         # number of fixed or free weights for 2nd-order latents

	# WHEN FORMATIVE INDICATORS, LOADINGS = 0
	for (j in 1:nlv1) {
		if (loadtype1[j] == 0) { C00[j,] = matrix(0,1,nvar) }	
	}
	
	C00 <- rbind(C00, matrix(0,nlv2,nvar))
	A00 <- cbind(C00, B00)

	W002 <- cbind(diag(1,nlv1), W002)

	V001 <- matrix(0,nvar,nvar)
	for (j in 1:nlv1) {
		nzaj <- which(C00[j,] != 0 )
		num_nzaj <- length(nzaj)
		for (i in 1:num_nzaj) {
			V001[nzaj[i],nzaj[i]] = 1
		}
	}
	
	# ---------------------------------------
	# model specification for all groups & initial random starts
	# ---------------------------------------
	
	W1i <- W001
	W2i <- W002
	Ai <- A00
	windex10 <- which(W001 == 99)
	windex20 <- which(W002 == 99)
	aindex0 <- which(A00 == 99)
	W01 <- matrix(0,ng*nvar,ng*nlv1)
	W02 <- matrix(0,ng*nlv1,ng*nlv)
	A0 <- matrix(0,ng*nlv,ntv)
	V0 <- matrix(0,ng*nvar,ng*ntv)
	W1 <- W01
	W2 <- W02
	A <- A0
	V <- V0
	kk <- 0
	ii <- 0
	ss <- 0
	tt <- 0
	for (j in 1:ng) {
		k <- kk + 1
		kk <- kk + nlv1
		i <- ii + 1
		ii <- ii + nvar
		s <- ss + 1
		ss <- ss + nlv
		t <- tt + 1
		tt <- tt + ntv
		W01[i:ii,k:kk] <- W001
		W1i[windex10] <- runif(length(windex10))
		W1[i:ii,k:kk] <- W1i
		W02[k:kk,s:ss] <- W002
		W2i[windex20] <- runif(length(windex20))
		W2[k:kk,s:ss] <- W2i
		A0[s:ss,] <- A00
		Ai[aindex0] <- runif(length(aindex0))
		A[s:ss,] <- Ai
		V0[i:ii,t:tt] <- cbind(V001,W001%*%W002)
		V[i:ii,t:tt] <- cbind(V001,W1i%*%W2i)
	}
	I <- matrix(0,ntv*ng,ntv)
	kk <- 0
	for (g in 1:ng) {
	    k <- kk + 1
		kk <- kk + ntv
		I[k:kk,] <- diag(1,ntv)
	}
	
	# generate orthogonal projector of equality constraints
	output.constmat <- constmat(A0)
	PHT <- output.constmat$PHT
	num_nzct <- output.constmat$num_nzct
	num_const <- output.constmat$num_const
	
	# ---------------------------------------
	# bootstrap starts here
	# ---------------------------------------
	num_par_W001 <- length(W001[!W001 == 0])
	num_nnz_C00 <- length(C00[!C00 == 0])
	num_nnz_B00 <- length(B00[!B00 == 0])

	vec_FIT <- matrix(0,nbt,1)
	vec_FIT_m <- matrix(0,nbt,1)
	vec_FIT_s <- matrix(0,nbt,1)
	vec_AFIT <- matrix(0,nbt,1)
	vec_GFI <- matrix(0,nbt,1)
	vec_SRMR <- matrix(0,nbt,1)

	MatW1 <- matrix(0,nbt,num_par_W001*ng)
	MatW2 <- matrix(0,nbt,num_par_W002*ng)
	Matload <- matrix(0,nbt,num_nnz_C00*ng)
	Matbeta <- matrix(0,nbt,num_nnz_B00*ng)
	Matsmc <- matrix(0,nbt,num_nnz_C00*ng)
	MatcorF <- matrix(0,nbt, nlv^2*ng)
	
	MatTE_S <- c()
	MatID_S <- c()
	MatTE_M <- c()
	MatID_M <- c()

	for (b in 0:nbt) {
	
		# generate a bootstrap sample Z (when b == 0, use the original data)
		if (b == 0) {
			if (moption > 1) {
				output.bootsample.imp <- bootsample.imp(z0, case_index, nvar, nobs_g, ng, b, nobs_tot, moption, missingvalue)
				Z <- output.bootsample.imp$Z
				z0_meanimp <- output.bootsample.imp$z0_meanimp
				rawz0 <- output.bootsample.imp$rawz0
			} else {
				output.bootsample <- bootsample(z0, case_index, nvar, nobs_g, ng, b, nobs_tot)
				Z <- output.bootsample$Z
			}
		} else {
			output.bootsample <- bootsample(z0, case_index, nvar, nobs_g, ng, b, nobs_tot)
			Z <- output.bootsample$Z
		}
	
		# parameter estimation
		if (b == 0) {
			if (moption == 3) {
				output.als.mg.ho2.imp <- als.mg.ho2.imp(Z, rawz0, W01, W02, A0, W1, W2, A, V, I, PHT, nvar, nlv1, nlv, ng, missingvalue, itmax, ceps)
				W <- output.als.mg.ho2.imp$W
				W1 <- output.als.mg.ho2.imp$W1
				W2 <- output.als.mg.ho2.imp$W2
				A <- output.als.mg.ho2.imp$A
				ZZ <- output.als.mg.ho2.imp$Z
				Psi <- output.als.mg.ho2.imp$Psi
				Gamma <- output.als.mg.ho2.imp$Gamma
				f <- output.als.mg.ho2.imp$f
				it <- output.als.mg.ho2.imp$it
				imp <- output.als.mg.ho2.imp$imp
			} else {
				output.als.mg.ho2 <- als.mg.ho2(Z, W01, W02, A0, W1, W2, A, V, I, PHT, nvar, nlv1, nlv, ng, itmax, ceps)
				W <- output.als.mg.ho2$W
				W1 <- output.als.mg.ho2$W1
				W2 <- output.als.mg.ho2$W2
				A <- output.als.mg.ho2$A
				Psi <- output.als.mg.ho2$Psi
				Gamma <- output.als.mg.ho2$Gamma
				f <- output.als.mg.ho2$f
				it <- output.als.mg.ho2$it
				imp <- output.als.mg.ho2$imp
			}
		} else {
			output.als.mg.ho2 <- als.mg.ho2(Z, W01, W02, A0, W1, W2, A, V, I, PHT, nvar, nlv1, nlv, ng, itmax, ceps)
			W <- output.als.mg.ho2$W
			W1 <- output.als.mg.ho2$W1
			W2 <- output.als.mg.ho2$W2
			A <- output.als.mg.ho2$A
			Psi <- output.als.mg.ho2$Psi
			Gamma <- output.als.mg.ho2$Gamma
			f <- output.als.mg.ho2$f
			it <- output.als.mg.ho2$it
			imp <- output.als.mg.ho2$imp
		}
		
		corF <- t(Gamma)%*%Gamma	# latent correlations across groups
		CR <- t(A[,1:nvar])			# loadings
		BR <- t(A[,(nvar+1):ntv])	# path coefficients

		DF <- nobs_tot*nvar
		npw1 <- length(which(W01 == 99))
		npw2 <- length(which(W02 == 99))
		dpht <- diag(PHT)
		if ( num_nzct == 0 ) {
			cnzt <- length(which(dpht == 1))
		} else {
			cnzt <- num_const + length(which(dpht == 1))
		}
		NPAR <- cnzt + npw1 + npw2	# number of parameters

		# model fit measures
		Fit <- 1 - f/sum(diag(t(Psi)%*%Psi))
		dif_m <- Psi[,1:nvar] - Gamma%*%t(CR)
		dif_s <- Psi[,(nvar+1):ntv] - Gamma%*%t(BR)
		Fit_m <- 1 - sum(diag(t(dif_m)%*%dif_m))/sum(diag(t(Z)%*%Z))
		Fit_s <- 1 - sum(diag(t(dif_s)%*%dif_s))/sum(diag(t(Gamma)%*%Gamma))
		Afit <- 1 - ((1-Fit)*(DF)/(DF - NPAR))
		output.modelfit.mg <- modelfit.mg(Z, W, A, nvar, nlv, ng, case_index)
		Gfi <- output.modelfit.mg$GFI
		Srmr <- output.modelfit.mg$SRMR
		COR_RES <- output.modelfit.mg$COR_RES

		# total and indirect effects in structural & measurement models
		total_s <- matrix(0,ng*nlv,nlv)			# total effects of latents in structural model
		indirect_s <- matrix(0,ng*nlv,nlv)		# indirect effects in structural model
		total_m <- matrix(0,ng*nlv,nvar)		# total effects of latent variables on indcators in measurement model
		indirect_m <- matrix(0,ng*nlv,nvar)		# indirect effect in measurement model
		k <- kk <- 0
		for (g in 1:ng) {
			k = kk + 1
			kk = kk + nlv
			output.effects <- effects(BR[,k:kk],CR[,k:kk])
			te_s <- output.effects$te_s
			ie_s <- output.effects$ie_s
			te_m <- output.effects$te_m
			ie_m <- output.effects$ie_m
			total_s[k:kk,] <- te_s
			indirect_s[k:kk,] <- ie_s
			total_m[k:kk,] <- te_m
			indirect_m[k:kk,] <- ie_m
		}

		# original sample solution
		if ( b == 0 ) {
		
			if ( moption == 2 ) { z0 <- z0_meanimp
			} else if ( moption == 3 ) {
				z0 <- matrix(0,nobs_tot,nvar)
				kk <- 0
				for (g in 1:ng) {
					k <- kk + 1
					kk <- kk + nvar
					z0[case_index[g,1]:case_index[g,2],] = ZZ[case_index[g,1]:case_index[g,2],k:kk]
					# ZZ - data after LS imputation
				}
			}
			
			if ( it <= itmax ) {
				if ( imp <= ceps ) {
					cat("The ALS algorithm converged in", it, "iterations (convergence criterion =", ceps, ")", "\n")
				} else {
					cat("The ALS algorithm failed to converge in", it, "iterations (convergence criterion =", ceps, ")", "\n")
				}
			}

			WR1 <- W1
			WR2 <- W2
			Cr <- CR
			Br <- BR
			samplesizes <- nobs_g	# sample size per group
			NPAR
			FIT <- Fit
			FIT_M <- Fit_m
			FIT_S <- Fit_s
			AFIT <- Afit
			GFI <- Gfi
			SRMR <- Srmr
			R2 <- matrix(0,ng,nlv)
			AVE <- matrix(0,ng,nlv1)
			Alpha <- matrix(0,ng,nlv1)
			rho <- matrix(0,ng,nlv1)      			# Dillong-Goldstein's rho (composite reliablity)
			Dimension <- matrix(0,ng,nlv1)			# Dimensionality per block
			lvmean <- matrix(0,ng,nlv)
			lvvar <- matrix(0,ng,nlv)
			corr_corres <- matrix(0,ng*nvar,nvar)	# matrix with correlations in lower trigular and correlation residuals in upper triangular
			ss <- 0
			kk <- 0
			for (g in 1:ng) {
				s <- ss + 1
				ss <- ss + nlv
				k <- kk + 1
				kk <- kk + nvar
				if ( moption == 3 ) {
					z0_g <- z0_meanimp[case_index[g,1]:case_index[g,2],]
				} else {
					z0_g <- z0[case_index[g,1]:case_index[g,2],]
				}
				W_g <- W[k:kk,s:ss]
				CF_g <- corF[s:ss,s:ss]
				B <- t(BR[,s:ss])
				for (j in 1:nlv) { R2[g,j] = t(B[,j])%*%CF_g[,j] }
				# quality measures of 1st-order latents for indicators  
				stdL <- CR[,s:ss]
				j2 <- 0
				for (j in 1:nlv1) {
					nnzload <- length(C00[j,][!C00[j,] == 0])
					j1 <- j2 + 1
					j2 <- j2 + nnzload
					if ( nnzload > 0 ) {
						sumload <- sum(stdL[,j]^2)
						sumload_rho1 <- sum(stdL[,j])^2
						sumload_rho2 <- sum(1-stdL[j1:j2,j]^2)
						AVE[g,j] <- sumload/nnzload
						rho[g,j] <- sumload_rho1/(sumload_rho1 + sumload_rho2)
					}
					zind = which(W001[,j] != 0)
					nzj = length(zind)
					if ( nzj > 1 ) {
						zsubset <- z0_g[,zind]
						Alpha[g,j] <- cronbach.alpha(zsubset)
						eigval <- svd(cor(zsubset))$d
						kr <- length(which(eigval>1))	# number of eigenvalues greater than 1
						Dimension[g,j] <- kr
					} else {
						Alpha[g,j] <- 1
					}
				}
				
				# calculate latent scores
				lvscore_g <- lvscore(z0_g, W_g)
				lvmean[g,] <- apply(lvscore_g,2,mean)
				lvvar[g,] <- apply(lvscore_g,2,var)
				sample_corr <- cor(z0_g)
				corr_corres[k:kk,][upper.tri(corr_corres[k:kk,], diag = FALSE)] <- COR_RES[k:kk,][upper.tri(COR_RES[k:kk,], diag = FALSE)]
				corr_corres[k:kk,][lower.tri(corr_corres[k:kk,], diag = FALSE)] <- sample_corr[lower.tri(sample_corr, diag = FALSE)]	
			}
			
			R2
			AVE
			Alpha
			rho
			LV_MEAN <- lvmean 	# Means of latent variables
			LV_VAR <- lvvar
			corr_corres
			Dimension
			mW1 <- as.matrix(W1[which(!W1 == 0)])
			mW2 <- as.matrix(W2[which(!W2 == 0)])
			mW2 <- as.matrix(mW2[mW2 != 1])
			mC <- as.matrix(CR[which(!CR == 0)])
			mB <- as.matrix(BR[which(!BR == 0)])
			mSMC <- mC^2
			mCF <- as.matrix(corF[which(!corF == 0)])
			latentcorr <- corF
			
			TE_S <- total_s
			ID_S <- indirect_s
			TE_M <- total_m
			ID_M <- indirect_m
	
		} else {	# bootstrap sample solution
			vecw1 <- as.matrix(W1[which(!W1 == 0)])
			bW2 <- as.matrix(W2[which(!W2 == 0)])
			vecw2 <- as.matrix(bW2[bW2 != 1])
			vecload <- as.matrix(CR[which(!CR == 0)])
			vecbeta <- as.matrix(BR[which(!BR == 0)])
			veccorF <- as.matrix(corF[which(!corF == 0)])
			
			vec_FIT[b] <- Fit
			vec_FIT_m[b] <- Fit_m
			vec_FIT_s[b] <- Fit_s
			vec_AFIT[b] <- Afit
			vec_GFI[b] <- Gfi
			vec_SRMR[b] <- Srmr
			
			MatW1[b,] <- t(vecw1)
			MatW2[b,] <- t(vecw2)
			Matload[b,] <- t(vecload)
			Matbeta[b,] <- t(vecbeta)
			Matsmc[b,] <- t(vecload^2)
			MatcorF[b,] <- t(veccorF)
			
			MatTE_S <- rbind(MatTE_S,total_s[which(!total_s == 0)])
			MatID_S <- rbind(MatID_S,indirect_s[which(!indirect_s == 0)])
			MatTE_M <- rbind(MatTE_M,total_m[which(!total_m == 0)])
			MatID_M <- rbind(MatID_M,indirect_m[which(!indirect_m == 0)])	
		}
	}
	
	# display bootstrap GSCA output
	if ( nbt > 0 ) {
	
		lb <- ceiling(nbt*0.025)
		ub <- ceiling(nbt*0.975)
		sortFIT <- sort(vec_FIT)
		sortFIT_m <- sort(vec_FIT_m)
		sortFIT_s <- sort(vec_FIT_s)
		sortAFIT <- sort(vec_AFIT)
		sortGFI <- sort(vec_GFI)
		sortSRMR <- sort(vec_SRMR)
		sortw1 <- apply(MatW1,2,sort)
		sortw2 <- apply(MatW2,2,sort)
		sortload <- apply(Matload,2,sort)
		sortbeta <- apply(Matbeta,2,sort)
		sortsmc <- apply(Matsmc,2,sort)
		sortcorF <- apply(MatcorF,2,sort)
		sortte_s <- apply(MatTE_S,2,sort)
		sortid_s <- apply(MatID_S,2,sort)
		sortte_m <- apply(MatTE_M,2,sort)
		sortid_m <- apply(MatID_M,2,sort)
		
		output.gsca.mg.ho2 <- list(WR1 = WR1, WR2 = WR2, CR = Cr, BR = Br, samplesizes = samplesizes, NPAR = NPAR,
								FIT = FIT, FIT_M = FIT_M, FIT_S = FIT_S, AFIT = AFIT, GFI = GFI, SRMR = SRMR,
								R2 = R2, AVE = AVE, Alpha = Alpha, rho = rho, LV_MEAN = LV_MEAN, LV_VAR = LV_VAR,
								corr_corres = corr_corres, Dimension = Dimension, latentcorr = latentcorr,
								TE_S = TE_S, ID_S = ID_S, TE_M = TE_M, ID_M = ID_M,
								mW1 = mW1, mW2 = mW2, mC = mC, mB = mB, mSMC = mSMC, mCF = mCF,
								lb = lb, ub = ub, vec_FIT = vec_FIT, vec_FIT_m = vec_FIT_m, vec_FIT_s = vec_FIT_s, vec_AFIT = vec_AFIT,
								vec_GFI = vec_GFI, vec_SRMR = vec_SRMR, MatW1 = MatW1, MatW2 = MatW2, Matload = Matload, Matbeta = Matbeta,
								Matsmc = Matsmc, MatcorF = MatcorF, MatTE_S = MatTE_S, MatID_S = MatID_S, MatTE_M = MatTE_M, MatID_M = MatID_M,
								sortFIT = sortFIT, sortFIT_m = sortFIT_m, sortFIT_s = sortFIT_s,
								sortAFIT = sortAFIT, sortGFI = sortGFI, sortSRMR = sortSRMR, sortw1 = sortw1, sortw2 = sortw2,
								sortload = sortload, sortbeta = sortbeta, sortsmc = sortsmc, sortcorF = sortcorF,
								sortte_s = sortte_s, sortid_s = sortid_s, sortte_m = sortte_m, sortid_m = sortid_m)
		output.gsca.mg.ho2

	} else {
	
		output.gsca.mg.ho2 <- list(WR1 = WR1, WR2 = WR2, CR = Cr, BR = Br, samplesizes = samplesizes, NPAR = NPAR,
								FIT = FIT, FIT_M = FIT_M, FIT_S = FIT_S, AFIT = AFIT, GFI = GFI, SRMR = SRMR,
								R2 = R2, AVE = AVE, Alpha = Alpha, rho = rho, LV_MEAN = LV_MEAN, LV_VAR = LV_VAR,
								corr_corres = corr_corres, Dimension = Dimension, latentcorr = latentcorr,
								TE_S = TE_S, ID_S = ID_S, TE_M = TE_M, ID_M = ID_M,
								mW1 = mW1, mW2 = mW2, mC = mC, mB = mB, mSMC = mSMC, mCF = mCF)
		output.gsca.mg.ho2	
	
	}

}