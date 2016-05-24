# unitroot.R 
# Copyright (C) 2014 by Matthew Clegg

# Generic functions related to testing for unit roots

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


ur_power <- function (ur_test, a0=0, a1=.95, trend=0, n=250, nrep=10000, p.value=0.05, ...) {
	# Uses simulation to estimate the power of the unit root test
	# function ur_test for an AR(1) process with a0, a1, and trend, where
	# the sample size is n.  The number of repititions performed
	# is given by nrep.  The function ur_test is presumed to return
	# a list or structure containing a field named p.value.  trend
    # may either be a scalar or a vector of length nrep.  In the latter 
    # case, each replication of the test uses a different value of trend.

    if (length(trend) == 1) trend <- rep(trend, nrep)
    if (length(trend) != nrep) stop("trend must have length 1 or nrep")
#	pvalues <- replicate(nrep, ur_test(rar1(n, a0, a1, t), ...)$p.value)
    pvalues <- sapply(1:nrep, function(i) ur_test(rar1(n, a0, a1, trend[i]), ...)$p.value)
	sum(pvalues <= p.value) / length(pvalues)
}

ur_power_table <- function (ur_test, nrep=1000, p.value=0.05,
	a1=c(0.995, 0.99, 0.98, 0.97, 0.96, 0.95),
    trend=0,
	n=c(100, 250, 500, 750, 1000, 1250),
    ...) {
	# Constructs a table of power estimates for realistic values of a1
	do.row <- function(nv) sapply(a1, function(a1v) ur_power(ur_test, 0, a1v, trend, nv, nrep, p.value, ...))
	pt <- do.call("rbind", lapply(n, function(nv) do.row(nv)))
	pt <- as.data.frame(pt)
	rownames(pt) <- n
	colnames(pt) <- a1
	pt
}

ur_specificity <- function (ur_test, rv_func, nrep=1000, p.value=c(0.01,0.02,0.05,0.1,0.2,0.5)) {
	# Uses simulation to estimate the accuracy of a unit root test.
	# Specificity is defined to be the probability that the null hypothesis
	# is accepted when the null hypothesis is true.  If the test is
	# calibrated properly, this will be identical the p-value reported
	# by the test.  Thus, the purpose of this function is to check that
	# p-values are being calculated accurately.  On input, ur_test is
	# a function for performing a unit root test, and rv_func is a function
	# that generates random variates that satisfy the null hypothesis.
	# For each p.value in the list, returns the percentage of trials where
	# the reported p-value was less than or equal to the specified value.
    #
    # Example: ur_specificity(pp.test, function() cumsum(rnorm(100)))
	
	replicates <- replicate(nrep, ur_test(rv_func())$p.value)
	results <- sapply(p.value, function(p) sum(replicates <= p)/length(replicates))
	names(results) <- p.value
	results
}

adf_power <- function (a0=0, a1=0.95, trend=0, n=250, nrep=10000, p.value=0.05, k=1) {
	# Uses simulation to estimate the power of the augmented dickey-fuller
	# test with k lags
	ur_power (function(X) adf.test(X,k=k), a0=a0, a1=a1, trend=trend, n=n, nrep=nrep, p.value=p.value)
}

adf_power_table <- function (nrep=1000, p.value=0.05,
	a1=c(0.995, 0.99, 0.98, 0.97, 0.96, 0.95),
    trend=0,
	n=c(250, 500, 750, 1000, 1250),
	k=1) {
	# Constructs a table of power estimates of the augmented dickey-fuller
	# test with k lags
	ur_power_table(function(X) adf.test(X,k=k), nrep=nrep, p.value=p.value, a1=a1, trend=trend, n=n)
}

