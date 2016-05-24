# FREGAT (c) 2016 Gulnara R. Svishcheva & Nadezhda M. Belonogova, ICG SB RAS

pval.famSKATO <- function(Z) {

	Z <- t(t(Z$Z) * (Z$w))
	m1 <- dim(Z)[2]
	Q.all <- c()
	for (rh in rhos) {
		CORB <- matrix(NA, m1, m1)
		if (rh < 1) {
			CORB05 <- chol(diag(m1) * (1 - rh) + rh)  # (correlation matrix for betas) ^ 0.5
		} else { CORB05 <- matrix(1, m1, m1) / sqrt(m1) }  # case of famBT (rh = 1)
		Z1 <- Z %*% t(CORB05)
		Q05 <- crossprod(Z1, SIG_res)  # (t(Z1) %*% SIG_res) 
		Q <- sum(Q05 * Q05)
		Q.all <- c(Q.all, Q)
	}
	geno <- P11CholInvCo %*% Z  # making centralized and independent genotypes # Cor
	Q.res <- NULL
	Q <- rbind(Q.all, Q.res)
	out <- SKAT_Optimal_Get_Pvalue(nullmod$total.var * Q / 2, geno / sqrt(2), rhos, method, acc, lim)
	out$p.value[1]

}

pval.famSKAT <- function(Z) {

	w <- Z$w
	Z <- Z$Z
	if (kernel == 'linear.weighted') {
		Q05 <- w * crossprod(Z, SIG_res)  # W %*% t(G) %*% Sigmai %*% res
		Q <- sum(Q05 * Q05)
		KerMat <- t(Z) %*% SPS %*% Z
		KKK <- t(KerMat * w) * w
	} else {
		K <- lskmTest.GetKernel(Z, kernel, w, dim(Z)[1], dim(Z)[2])
		Q <- sum(SIG_res * crossprod(K, SIG_res))
		KKK <- CholSigmaiP11 %*% K %*% t(CholSigmaiP11)
	}

	eig <- eigen(KKK, symmetric = TRUE, only.values = TRUE)
	ev <- eig$values[eig$values > 1e-6 * eig$values[1]]

	if (method == 'kuonen') {
		p <- pchisqsum(Q, rep(1, length(ev)), ev, lower.tail = F, method = 'sad') }
	else if (method == 'davies') {
		p <- davies(Q, ev, acc = acc, lim = lim)$Qq
	}
	p <- max(min(p, 1), 0)

	if (!return.variance.explained) { p }

	else {

		if (!reml) {

			ZZ <- t(t(Z) - colMeans(Z))
			ZZ <- t(t(ZZ)* w)
			Ker <- ZZ %*% t(ZZ)
			param <- c(nullmod$total.var, nullmod$h2, 0.5, nullmod$alpha[match(colnames(X0)[-1], rownames(nullmod$alpha)), 1])
			fullik <- optim(param, tau3, X = X0[, -1], R = kin * 2, Ker = Ker, y = X0[, 1], I = diag(n1),
				lower = c(rep(0 + .Machine$double.eps, 3), rep(-Inf, length(param) - 3)),
				upper = c(Inf, rep(1 - .Machine$double.eps, 2), rep(Inf, length(param) - 3)), method ='L-BFGS-B')
			h2r <- fullik$par[2] * (1. - fullik$par[3]) # hk2 * (1. - dd) - proportion of variance by kernel (K)
			h2 <- fullik$par[2] * fullik$par[3] # hk2 * dd  - proportion of variance by relationship (R)
			if (h2r < 0) h2r <- 0
			if (h2 < 0) {
				h2 <- 0
				h2r <- NA
			}
			if (fullik$par[1] < 0) {
				fullik$par[1] <- 0
				h2 <- NA
				h2r <- NA
			}
			#fullik$value = (-2*log(LH)-n*log(2*pi))
			LH <- -(fullik$value + n * log(2 * pi)) # 2*logLH1
			chi <- LH - 2 * nullmod$logLH # 2*logLH1 - 2*logLH0
			p <- pchisq(chi, 1, lower.tail = F)
#				chi <- -fullik$value - n*log(2*pi) - 2*nullmod$logLH #   chi = 2*logLH1 - 2*logLH0
#	p <- pchisq(chi, 1, lower.tail = F)
			LH <- LH / 2 # logLH1
			#LH<- (-fullik$value - n*log(2*pi))/2
			c(p, h2r, h2, fullik$par[1], LH)

		} else { # reml estimate
			Z <- t(t(Z) * w)
			geno <- P11CholInvCo %*% Z # Cov
			A <- t(geno) %*% geno
			T <- eigen(A, symmetric = TRUE)  # T %*% lam %*% t(T)
			lam <- T$val
			zzz <- t(T$vec) %*% t(geno) %*% pheno  # vector
			a.est <- optimize(tau.fun, interval = c(0, 1), lam = lam, zzz = zzz)$minimum
			tau2 <- a.est/(1-a.est)
			eZ <- colMeans(Z)
			cf <- tau2 * (sum(Z * Z)- sum(eZ*eZ)*n1)
			prop <- cf / (sum(diag(SIGMA)) + cf)
			c(p, prop)#, tau2, cf)
		}
	}

}
