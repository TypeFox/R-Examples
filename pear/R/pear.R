`pear` <-
function(z, m, ic = "none")
{
	if(ic == "none") {
		if(missing(m)) {
			stop("Error: m required. Vector of model orders")
		}
		lag.max <- max(m)
		acf.out <- peacf(z, lag.max, plot = FALSE)
	}
	else if(!(ic == "none")) {
		p <- attr(z, "tsp")[3]
		N <- length(z)
		lag.max <- ceiling((0.25 * N)/p)
		pacf.out <- pepacf(z, lag.max, plot = FALSE)
		if(ic == "aic") {
			m <- pacf.out$maice
		}
		else if(ic == "bic") {
			m <- pacf.out$mbice
		}
		lag.max <- max(m)
		acf.out <- pacf.out$acf.out
	}
	p <- acf.out$period
	if((ic == "none") && (length(m) != p))
		m <- rep(m, length = p)
	res <- numeric(length(z))
	means <- acf.out$means
	for(imonth in 1:p) {
		k <- cycle(z) == imonth
		res[k] <- means[imonth]
	}
	res <- z - res
	attr(res, "tsp") <- attr(z, "tsp")
	z <- res
	Qm <- 0
	Qm.sl <- 0
	nyrs <- acf.out$sub.lengths
	bsd <- acf.out$benchmark.sd
	phi <- matrix(numeric(1), nrow = p, ncol = lag.max)
	se.phi <- matrix(numeric(1), nrow = p, ncol = lag.max)
	resvar <- numeric(p)
	sd <- acf.out$standard.deviations
	acvf <- cbind(rep(1, p), acf.out$acf)
	cov <- vector("list", p)
	for(i in 1:p) {
		for(j in 1:(lag.max + 1)) {
			acvf[i, j] <- acvf[i, j] * sd[i] * sd[((i - j) %% p) + 
				1]
		}
	}
	for(imonth in 1:p) {
		pm <- m[imonth]
		if(pm == 0) {
			resvar[imonth] <- acvf[imonth, 1]
			phim <- 0
			sephim <- 0
		}
		else {
			a <- matrix(numeric(1), nrow = pm, ncol = pm)
			for(i in 1:pm) {
				for(j in 1:pm) {
				  k <- i - j
				  kmonth <- (imonth - j - 1) %% p + 1
				  if(k < 0) {
				    k <-  - k
				    kmonth <- (imonth - i - 1) %% p + 1
				  }
				  a[i, j] <- acvf[kmonth, k + 1]
				}
			}
			b <- acvf[imonth, 1 + (1:pm)]
			phim <- solve(a, b)
			resvar[imonth] <- acvf[imonth, 1] - phim %*% acvf[
				imonth, 1 + (1:pm)]
			cov[[imonth]] <- (solve(a) * resvar[imonth])/nyrs[
				imonth]
			sephim <- sqrt(diag(cov[[imonth]]))
			if(pm < lag.max) {
				phim <- c(phim, rep(0, lag.max - pm))
				sephim <- c(sephim, rep(0, lag.max - pm))
			}
		}
		phi[imonth,  ] <- phim
		se.phi[imonth,  ] <- sephim
	}
	for(i in 1:length(z)) {
		res[i] <- 0
		imonth <- cycle(z)[i]
		pm <- m[imonth]
		if(((i - pm) > 0) && (pm > 0))
			res[i] <- z[i] - z[i - (1:pm)] %*% phi[imonth, 1:pm]
		else if(pm == 0)
			res[i] <- z[i]
	}
	ra <- peacf(res, plot = FALSE)
	residual.acf <- ra$acf
	dimnames(residual.acf)[[1]] <- dimnames(acf.out$acf)[[1]]
	L.racf <- ncol(residual.acf)
	residual.acf.sd <- matrix(numeric(1), nrow = p, ncol = L.racf, dimnames
		 = dimnames(residual.acf))
	psi <- cbind(rep(1, p), pepsi(phi, L.racf - 1))
	for(imonth in 1:p) {
		pm <- m[imonth]
		if(pm > 0) {
			X <- matrix(numeric(1), nrow = L.racf, ncol = pm)
			for(i in 1:L.racf) {
				for(j in 1:pm) {
				  if((i - j) >= 0) {
				    imonthmi <- (imonth - i - 1) %% p + 1
				    jmonth <- (imonth - j - 1) %% p + 1
				    rsigmas <- sqrt(resvar[imonthmi]/resvar[
				      imonth])
				    X[i, j] <- psi[jmonth, i - j + 1] * rsigmas
				  }
				}
			}
			covest <- cov[[imonth]]
			residual.acf.sd[imonth,  ] <- sqrt(abs(1/nyrs[imonth] - 
				diag(X %*% covest %*% t(X))))
		}
		else if(pm == 0) {
			residual.acf.sd[imonth,  ] <- sqrt(1/nyrs[imonth])
		}
	}
	dimnames(phi) <- list(dimnames(acf.out$acf)[[1]], paste("lag", 1:ncol(
		phi)))
	dimnames(se.phi) <- dimnames(phi)
	names(m) <- dimnames(phi)[[1]]
	names(resvar) <- names(m)
	names(cov) <- names(m)
	QM <- ra$portmanteau.test$QM
	QM.df <- ra$portmanteau.test$QM.df - matrix(m, nrow = p, ncol = ncol(QM
		))
	QM.sl <- matrix((1 - pchisq(QM, QM.df)), nrow = nrow(QM), ncol = ncol(
		QM), dimnames = dimnames(QM))
	portmanteau.test <- list(QM = QM, QM.df = QM.df, QM.sl = QM.sl)
	out <- list(model.orders = m, phi = phi, se.phi = se.phi, resvar = 
		resvar, residuals = res, portmanteau.test = portmanteau.test, 
		residual.acf = residual.acf, residual.acf.sd = residual.acf.sd, 
		cov = cov)
	out
}

