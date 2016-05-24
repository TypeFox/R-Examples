# egcm_base.R
# Copyright (C) 2014 by Matthew Clegg

# A collection of basic functions that support the egcm package.

#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/


YMD <- function (days=0) {
  # Returns today's date + days formatted as an integer in YYYYMMDD format.
  as.numeric(format(Sys.Date()+days, "%Y%m%d"));
}

demean <- function (Y, ...) {
    # Centers Y around its mean.
    Y - mean(Y, ...)
}

detrend <- function (Y) {
	# Removes a linear trend from Y
	if (is(Y, "zoo")) {
		X <- as.numeric(index(Y))
		X <- X - X[1]
		Y <- coredata(Y)
	} else {
		X <- 1:length(Y)
	}
	beta = cov(X,Y)/var(X)
	alpha = mean(Y) - beta * mean(X)
	eps = Y - alpha - beta * X
	eps
}

quantile_table_interpolate <- function (qtab, sample_size, stat, stop.on.na=FALSE) {
	# On input, qtab is a dataframe of quantiles.  Each column corresponds to
	# a sample size, and each row corresponds to a quantile value.  The sample
	# sizes are given in the first row, and the quantiles are given in the
	# first column.  
	n <- nrow(qtab)
	i <- findInterval(sample_size, qtab[1,2:ncol(qtab)])+1
	if (i == 1) {
		parent_name <- as.character(sys.call(-1)[[1]])
		if (stop.on.na) {
			stop (parent_name," requires a minimum of ", qtab[1,2], " observations.")
		} else {
			warning (parent_name, " requires a minimum of ", qtab[1,2], " observations.")
			return(NA)
		}
	}
	y1 <- approx(qtab[2:n, i], qtab[2:n, 1], stat, rule=2)$y
	if (i < ncol(qtab)) {
		y2 <- approx(qtab[2:n, i+1], qtab[2:n, 1], stat, rule=2)$y
		n1 <- qtab[1,i]
		n2 <- qtab[1,i+1]
		y <- y1 * (n2 - sample_size) / (n2 - n1) + y2 * (sample_size - n1)/(n2 - n1)
	} else {
		y <- y1
	}
	y
}

acor <- function (X, k=1, na.rm=FALSE) {
	# Calculates the lag k autocorrelation of X, e.g.,
	#   cov(X[t], X[t+k])/var(X)
	Xc <- coredata(X)
	n <- length(Xc)
	I1 <- 1:(n-k)
	I2 <- (k+1):n
	if (na.rm) {
		ac <- cov(Xc[I1], Xc[I2], use="complete.obs")/var(Xc[I1], na.rm=TRUE)
	} else {
		ac <- cov(Xc[I1], Xc[I2])/var(Xc[I1])
	}
	ac
}

###################################################
##### Functions for Generating Random Variates ####
###################################################

rar1 <- function (n, a0=0, a1=1, trend=0, sd=1, x0=0) {
	# Generates a vector of length n representing a simulation of an AR(1)
	# process   X[k] = a0 +  a1 * X[k-1] + eps
	# where eps is an i.i.d. normal random variate with mean 0 and standard
	# deviation sd.  
    #
    # If trend is non-zero, then returns a realization of a trend-stationary 
    # AR(1) process.  E.g., the process is defined by the relations:
    #    R[k] = a0 + a1 * R[k-1] + eps
    #    X[k] = k * trend + R[k]
	eps <- rnorm(n, 0, sd)
	X <- numeric(n)
	xp <- x0
	for (k in 1:n) {
		X[k] <- xp <- a0 + a1 * xp + eps[k]
	}
	X + trend * (1:n)
}

rcoint <- function (n, alpha=runif(1,-10,10), beta=runif(1,-10,10), rho=runif(1,0,1),
  sd_eps = 1, sd_delta=1, X0 = 0, Y0 = 0) {
	# Generates a random pair of cointegrated vectors X[t] and Y[t] subject 
	# to the relations:
	#   X[t] = X[t-1] + eps[t]
	#   Y[t] = alpha + beta * X[t] + R[t]
	#   R[t] = a2 * R[t-1] + delta[t]
	# where eps, delta are NID(0,1).  Returns the n x 2 matrix containing 
	# X and Y.
	X <- rep(0, n+1)
	Y <- rep(0, n+1)
	eps <- rnorm(n, 0, sd_eps)
	delta <- rnorm(n, 0, sd_delta)
    X[1] <- X0
    Y[1] <- Y0
	R <- Y0 - alpha - beta * X0
	for (i in 2:(n+1)) {
		X[i] <- X[i-1] + eps[i-1]
		R <- rho * R + delta[i-1]
		Y[i] <- alpha + beta * X[i] + R
	}
	M <- cbind(X[2:(n+1)],Y[2:(n+1)])
	colnames(M) <- c("X", "Y")
	M
}

