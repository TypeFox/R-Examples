`pepacf` <-
function(z, lag.max, plot = TRUE, acf.out)
{
	if(missing(acf.out)) {
		if(missing(lag.max))
			acf.out <- peacf(z, plot = FALSE)
		else acf.out <- peacf(z, lag.max, plot = FALSE)
	}
	p <- acf.out$period
	if(missing(lag.max))
		lag.max <- ceiling(length(z)/(p * 4))
	nyrs <- acf.out$sub.lengths
	bsd <- acf.out$benchmark.sd
	w1 <- numeric(p)
	w2 <- numeric(p)
	w3 <- numeric(p)
	pacf <- matrix(numeric(1), nrow = p, ncol = lag.max)
	phi <- matrix(numeric(1), nrow = p, ncol = lag.max)
	phiF <- matrix(numeric(1), nrow = p, ncol = lag.max)
	resvar <- matrix(numeric(1), nrow = p, ncol = lag.max + 1)
	resvarF <- matrix(numeric(1), nrow = p, ncol = lag.max + 1)
	sd <- acf.out$standard.deviations
	acvf <- cbind(rep(1, p), acf.out$acf)
	for(i in 1:p) {
		for(j in 1:(lag.max + 1)) {
			acvf[i, j] <- acvf[i, j] * sd[i] * sd[((i - j) %% p) + 
				1]
		}
	}
	resvar[, 1] <- acvf[, 1]
	resvarF[, 1] <- acvf[, 1]
	for(ilag in 1:lag.max) {
		SFPREV <- resvarF[p, ilag]
		for(imonth in 1:p) {
			D <- acvf[imonth, ilag + 1]
			if(ilag > 1) {
				for(i in 1:(ilag - 1)) {
				  mmi <- ((imonth - i - 1) %% p) + 1
				  D <- D - acvf[mmi, (ilag - i + 1)] * phi[
				    imonth, i]
				}
			}
			phi[imonth, ilag] <- D/SFPREV
			phiF[imonth, ilag] <- D/resvar[imonth, ilag]
			pacf[imonth, ilag] <- D/sqrt(SFPREV * resvar[imonth, 
				ilag])
			TEMP <- 1 - phi[imonth, ilag] * phiF[imonth, ilag]
			resvar[imonth, ilag + 1] <- resvar[imonth, ilag] * TEMP
			SFKEEP <- resvarF[imonth, ilag]
			resvarF[imonth, ilag + 1] <- SFPREV * TEMP
			SFPREV <- SFKEEP
			if(ilag > 1) {
				for(i in 1:(ilag - 1)) {
				  w2[i] <- w1[i] - phi[imonth, ilag - i] * phiF[
				    imonth, ilag]
				  w3[i] <- phi[imonth, i] - phi[imonth, ilag] * 
				    w1[ilag - i]
				}
				for(i in 1:(ilag - 1)) {
				  w1[i] <- phiF[imonth, i]
				  phiF[imonth, i] <- w2[i]
				  phi[imonth, i] <- w3[i]
				}
			}
		}
		w1 <- phiF[p,  ]
	}
	nyrs.matrix <- matrix(nyrs, ncol = lag.max + 1, nrow = p)
	aic <- nyrs.matrix * log(resvar)
	bic <- aic
	k.matrix <- matrix(0:lag.max, byrow = TRUE, nrow = p, ncol = lag.max + 1)
	aic <- aic + 2 * k.matrix
	bic <- bic + log(nyrs.matrix) * k.matrix
	maice <- find.ice(aic)
	mbice <- find.ice(bic)
	r <- list(acf.out = acf.out, pacf = pacf, residual.sd = sqrt(resvar), 
		phi = phi, aic = aic, bic = bic, maice = maice, mbice = mbice)
	attr(r, "type") <- "pacf"
	if(plot) {
		peacf.plot(r)
	}
	r
}

