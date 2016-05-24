`pepsi` <-
function(phi, lag.max)
{
	phi[is.na(phi)] <- 0
	p <- nrow(phi)
	m.max <- ncol(phi)
	if(m.max > 1)
		m <- t(apply(t(apply(t(apply(phi, 1, "rev")) != 0, 1, "cumsum")
			) != 0, 1, "cumsum"))[, ncol(phi)]
	else m <- rep(1, p)
	psi <- matrix(numeric(1), nrow = p, ncol = lag.max + 1)
	psi[, 1] <- 1
	for(ilag in 1:lag.max) {
		for(imonth in 1:p) {
			TEMP <- numeric(1)
			for(jlag in 1:m[imonth]) {
				kmonth <- (imonth - jlag - 1) %% p + 1
				klag <- ilag - jlag
				if(klag >= 0)
				  TEMP <- TEMP + phi[imonth, jlag] * psi[kmonth,
				    klag + 1]
			}
			psi[imonth, ilag + 1] <- TEMP
		}
	}
	psi <- psi[, -1]
	dimnames(psi) <- list(periods = paste("period", 1:p), lags = paste(
		"lag", 1:lag.max))
	psi
}

