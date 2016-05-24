`peacf` <-
function(z, lag.max, plot = TRUE)
{
	if(!is.ts(z)) {
		stop("Error: input is not a time series.\nUse ts() to fix-up.")
	}
	p <- attr(z, "tsp")[3]
	if(missing(lag.max))
		lag.max <- ceiling(length(z)/(p * 4))
	z.names <- attr(z, "period.abb")
	zc <- z
	zmean <- numeric(p)
	zsd <- numeric(p)
	nyrs <- numeric(p)
	nyrs.max <- ceiling(length(z)/p)
	for(imonth in 1:p) {
		k <- cycle(z) == imonth
		zmean[imonth] <- mean(z[k])
		zc[k] <- zmean[imonth]
		nyrs[imonth] <- length(z[k])
		zsd[imonth] <- sqrt(sum((z[k] - zc[k])^2)/nyrs.max)
	}
	zc <- z - zc
	if(lag.max > min(nyrs)) {
		cat("\nwarning: lag.max exceeds data span so reset to a smaller value"
			)
		lag.max <- min(nyrs) - 1
	}
	periodic.acf <- matrix(numeric(1), nrow = p, ncol = lag.max)
	for(ilag in 1:lag.max) {
		zac1 <- ts.intersect(zc, lag(zc,  - ilag))
		zac <- ts(apply(zac1, MARGIN = 1, FUN = prod), frequency = p, 
			start = attr(zac1, "tsp")[1])
		for(imonth in 1:p) {
			imonthlag <- (imonth - ilag - 1) %% p + 1
			periodic.acf[imonth, ilag] <- sum(zac[cycle(zac) == 
				imonth])/(nyrs.max * zsd[imonth] * zsd[
				imonthlag])
		}
	}
	if(is.null(z.names))
		dimnames(periodic.acf) <- list(periods = paste("period", 1:p), 
			lags = paste("lag", 1:lag.max))
	else dimnames(periodic.acf) <- list(periods = z.names, lags = paste(
			"lag", 1:lag.max))
	bsd <- 1/sqrt(max(nyrs))
	title.z <- attr(z, "title")
	if(is.null(title.z)) {
		title.z <- " "
	}
	Q1 <- nyrs.max * sum(periodic.acf[, 1]^2)
	Q1.sl <- 1 - pchisq(Q1, p)
	periodicity.test <- list(Q1 = Q1, Q1.sl = Q1.sl)
	null.var <- matrix(numeric(1), nrow = p, ncol = lag.max)
	for(i in 1:p) {
		for(j in 1:lag.max) {
			null.var[i, j] <- var.periodic.correlation(j, i, nyrs[i
				], p)
		}
	}
	QM <- matrix(t(apply((periodic.acf^2)/null.var, 1, cumsum)), nrow = p, 
		ncol = ncol(periodic.acf))
	QM.df <- col(QM)
	if(lag.max >= 5)
		QM.lags <- if(0 == (lag.max %% 5)) 5 * (1:(lag.max/5)) else c(5 *
				(1:trunc(lag.max/5)), lag.max)
	else QM.lags <- lag.max
	QM <- matrix(QM[, QM.lags], nrow = p, ncol = length(QM.lags))
	QM.df <- matrix(QM.df[, QM.lags], nrow = p, ncol = length(QM.lags))
	portmanteau.test <- list(QM = QM, QM.df = QM.df)
	r <- list(means = zmean, standard.deviations = zsd, acf = periodic.acf, 
		benchmark.sd = bsd, sub.lengths = nyrs, period = p, title = 
		title.z, periodicity.test = periodicity.test, portmanteau.test
		 = portmanteau.test)
	attr(r, "type") <- "acf"
	if(plot) {
		peacf.plot(r)
	}
	r
}

