# egcm.R  
# Copyright (C) 2014 by Matthew Clegg

#  Engle-Granger Cointegration Models for Pairs Trading

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

# The purpose of this module is to provide a simple implementation of 
# the Engle Granger cointegration model that is convenient for use in
# the analysis of pairs trades.  
#
# Given two time series Y[t] and X[t], the Engle Granger cointegration 
# model in its simplest form looks for alpha, beta and rho such that
#
#    Y[t] = alpha + beta * X[t] + R[t]
#    R[t] = rho * R[t-1] + epsilon[t]
#
# where epsilon[t] is a series of independent and identically distributed
# innovations with mean zero.  If alpha, beta and rho can be found such that
# -1 < rho < 1, then the series are said to be cointegrated.  If abs(rho) = 1,
# then the residual series R[t] is said to have a unit root (or alternatively, 
# it is said to follow a random walk).
#
# It should be noted that Engle and Granger's 1987 paper allows for a more general
# definition of cointegration.  Namely, there can be multiple series included
# in the cointegrating relationship.  Also, the residual series R[t] may
# be any ARMA process.  If this greater generality is needed, the reader
# is referred to the urca package of Bernhard Pfaff.
#
# The major functions provided by this module are as follows:
#
# egcm(X,Y)          -- Constructs an Engle-Granger cointegration model from X & Y
# summary.egcm(E)    -- Prints various summary statistics on the Engle-Granger
#                       cointegration model constructed from X & Y
# plot.egcm(E)       -- Creates a graph of the Engle-Granger cointegration model
#
# The following ancillary functions are also provided:
#
# rcoint(n, alpha, beta, rho) -- Generates a random pair of cointegrated vectors
# egc_test_specificity()   -- Calculates the specificity of a cointegration (unit root) test
# egc_test_power()         -- Calculates the power of a cointegration (unit root) test
# egc_test_power_table()   -- Calculates a table of powers of a cointegration (unit root) test
# as.data.frame.egcm()     -- Converts an egcm object to a single row data.frame
# residuals.egcm()         -- Returns the residual series R[t] associated to an egcm object
# innovations.egcm()       -- Returns the series epsilon[t] of innovations associated to an
#                             egcm object
# test.egcm.innovations()  -- Tests the goodness of fit of a set of innovations
# 
#
# References
#
# Engle, R. F., & Granger, C. W. (1987). Co-integration and error correction: 
#   representation, estimation, and testing. Econometrica, (55) 2, 251-276.
# Pfaff, Bernhard (2008).  Analysis of Integrated and Cointegrated Time 
#   Series with R.  Springer.


################################################
##### Engle Granger Cointegration Test     #####
################################################

# The EGCM package maintains several writable variables that are local to the
# package.  These are saved in the following environment and accessed with getter
# and setter methods.
egcm.env <- new.env ()

# There are so many different unit root tests that have been implemented that it
# seems easier to create a global variable enumerating them.

egcm.urtests.internal <- c(
	"adfraw",    # Augmented Dickey-Fuller (adf.test{tseries}), using default p-values
	
# All of the tests below use p-values that have been re-calibrated
	"adf",       # Augmented Dickey-Fuller test (adf.test{tseries})
	"pp",        # Phillips-Perron test (pp.test{tseries})
	"jo-e",      # Johansen's eigenvalue test (ca.jo{urca})
	"jo-t",      # Johansen's trace test (ca.jo{urca})
	"ers-p",	 # Elliott, Rothenberg and Stock point optimal test (ur.ers{urca})
	"ers-d",     # Elliott, Rothenberg and Stock DF-GLS test (ur.ers{urca})
	"sp-r",      # Schmidt and Phillips rho statistics (ur.sp{urca})
	"hurst",     # Hurst exponent (aggvarFit{fArma})
	"bvr",       # Breitung's variance ratio
	"pgff"       # Pantula, Gonzales-Farias and Fuller Rho statistic 
)

egcm.set.default.urtest <- function (urtest) {
	# Sets the default unit root test that will be used in subsequent calls
	# to test.egcm and egcm.
	
	if (!(urtest %in% egcm.urtests.internal)) {
		stop ("No such unit root test: ", urtest)
	}
    
#    if (urtest %in% c("jo-e", "jo-t", "ers-p", "ers-d", "sp-r")) require(urca)
#    if (urtest == "hurst") require(fArma)

	urtests.list <<- c(urtest, setdiff(egcm.urtests.internal, urtest))
    assign("urtest.default", urtests.list, envir=egcm.env)
}

egcm.urtests <- function () {
    # Returns a list of the mnemonics of all available unit root tests
    get("urtest.default", envir=egcm.env)
}

egcm.default.urtest <- function() {
    # Returns the default unit root test that is used by test.egcm and egcm
    egcm.urtests()[1]
}

# egcm.set.default.urtest("pp")

# The following global variable enumerates the tests that are available for
# checking if the input series are integrated.  The tests from the urca library
# have been omitted because these tests do not provide p-values.

egcm.i1tests.internal <- c(
	"adf",	 	# Augmented Dickey-Fuller (adf.test{tseries}), default p.values
	"pp",       # Phillips-Perron test (pp.test{tseries}), default p.values
	"bvr",      # Breitung's variance ratio
	"pgff"      # Pantula, Gonzales-Farias and Fuller Rho statistic
)

egcm.set.default.i1test <- function (i1test) {
	# Sets the default test for integration that will be used in subsequent calls to egcm.
	
	if (!(i1test %in% egcm.i1tests.internal)) {
		stop ("No such integration test: ", i1test)
	}

	i1tests.list <- c(i1test, setdiff(egcm.i1tests.internal, i1test))
    assign("i1test.default", i1tests.list, envir=egcm.env)
}

egcm.i1tests <- function () {
    # Returns a list of all of the mnemonics of all available I(1) tests,
    # with the currently selected default at the head of the list.
    get("i1test.default", envir=egcm.env)
}

egcm.default.i1test <- function () {
    # Returns the default I(1) test that is used by egcm.
    egcm.i1tests()[1]
}

# egcm.set.default.i1test("pp")

# The following getter/setter functions are used to control the
# confidence level (p-value) that is used for the various statistical
# tests.

egcm.set.default.pvalue <- function (p) {
    # Sets the default p-value used in the tests of cointegration
    if (p < 0.001) stop("p values less than 0.001 are not supported")
    assign("pvalue.default", p, envir=egcm.env)
}

egcm.default.pvalue <- function () {
    # Returns the default p-value used in the tests of cointegration
    get("pvalue.default", envir=egcm.env)
}

# Set the following to TRUE if you would like egcm.summary to print
# its estimates of the expected investment returns that would be
# obtained from a pairs trade.  This is an experimental feature that
# should not be relied upon in practice.  USE AT YOUR OWN RISK!

include.investment.scenarios <- function (v) {
    if (missing(v)) {
        get("include.investment.scenarios", envir=egcm.env)
    } else {
        if(!is.logical(v)) stop("argument to include.investment.scenarios must be logical. ",
            "Instead received ", v)
        assign("include.investment.scenarios", v, envir=egcm.env)
    }
}

include.investment.scenarios(FALSE)

test.egcm <- function(EGCM, test.method=egcm.default.urtest()) {
	# Tests whether EGCM is cointegrated by performing a unit root test on the
	# residual series.  The choice of unit root test is determined by test.method:
	#     adfraw: Augmented Dickey-Fuller test as implemented by adf.test{tseries}
	#     adf:   Augmented Dickey-Fuller test, as implemeneted by adf.test{tseries}.
	#     pp:    Phillips-Perron test, as implemented by pp.test{tseries}.
	#     pgff:  Weighted symmetric estimator rho of Pantula, Gonzales-Farias and Fuller.
	#     jo-e:  Eigenvalue test from Johansen's VECM model as implemented by ca.jo{urca}.
	#     jo-t:  Trace test from Johansen's VECM model as implemented by ca.jo{urca}.
	#     ers-p: Point optimal test of Elliott, Rothenberg and Stock as implemented by ur.ers{urca}.
	#     ers-d: DF-GLS test of Elliott, Rothenberg and Stock as implemented by ur.ers{urca}.
	#     sp-r:  Rho statistic from Schmidt and Phillips, as implemented by ur.sp{urca}.
	#     bvr:   Breitung's variance ratio test.
	#     hurst: Hurst exponent as calculated by aggvarFit{fArma}.
	#
	# EGCM can either be an egcm object returned by egcm(), or it can be a two-column
	# matrix.  In the latter case, a regression is performed of the first column on the
	# second to obtain the residual series, which is tested for a unit root.
	#
	# Returns an object of type htest representing the results of the hypothesis test.
	# A low p.value is interpreted as evidence that the null hypothesis
	# of a unit root should be rejected.  In other words, a low p.value is evidence
    # for cointegration.

	test.method <- match.arg(test.method, egcm.urtests())
    DNAME <- deparse(substitute(EGCM))

#    if (test.method %in% c("jo-e", "jo-t", "ers-p", "ers-d", "sp-r")) require(urca)
#    if (test.method == "hurst") require(fArma)
	
	if (is(EGCM, "egcm")) {
		R <- EGCM$residuals
		X <- cbind(EGCM$S1, EGCM$S2)
	} else {
		X <- coredata(EGCM)
		beta <- cov(X[,2], X[,1]) / var(X[,1])
		alpha <- mean(X[,2]) - beta * mean(X[,1])
		R <- X[,2] - beta * X[,1] - alpha		
	}
	
	if (test.method %in% c("jo-e", "jo-t")) {
		if (test.method == "jo-e") {
			jo <- ca.jo(provideDimnames(X), type="eigen", ecdet="const")
			STAT <- -jo@teststat[2]
			names(STAT) <- "lambda_max"
			PVAL <- quantile_table_interpolate(egc_joe_qtab, length(R), STAT)
			URTEST <- jo@test.name		
		} else if (test.method == "jo-t") {
			jo <- ca.jo(provideDimnames(X), type="trace", ecdet="const")
			STAT <- -jo@teststat[2]
			names(STAT) <- "trace"
			PVAL <- quantile_table_interpolate(egc_jot_qtab, length(R), STAT)
			URTEST <- jo@test.name		
		}
	    htest <- structure(list(statistic = STAT, alternative = "cointegrated", 
	        p.value = PVAL, method = "", urtest = URTEST, data.name = DNAME), 
	        class = "htest")    	    
	} else {
		htest <- egc.residuals.test(R, test.method)
		htest$data.name <- DNAME
	}

	htest$method <- sprintf("Engle Granger two-step cointegration test (%s)", test.method)	
	htest
}

egc.residuals.test <- function(R, test.method=setdiff(egcm.urtests(), c("jo-e", "jo-t"))) {
	# Tests whether the residual series from the Engle-Granger procedure
	# contains a unit root.  
	#
	# Input values:
	#   R:    The residual series that is to be tested
	#   test.method:  The method to be used for testing for a unit root.
	#     One of the following choices is permitted:
	#
    #     adfraw: Augmented Dickey-Fuller test, as implemeneted by adf.test{tseries}.
	#     adf:   Augmented Dickey-Fuller test, but with re-calibrated p-values.
	#     pp:    Phillips-Perron test, as implemented by pp.test{tseries}.
	#     pgff:  Weighted symmetric estimator rho of Pantula, Gonzales-Farias and Fuller.
	#     ers-p: Point optimal test of Elliott, Rothenberg and Stock as implemented by ur.ers{urca}.
	#     ers-d: DF-GLS test of Elliott, Rothenberg and Stock as implemented by ur.ers{urca}.
	#     sp-r:  Rho statistic from Schmidt and Phillips, as implemented by ur.sp{urca}.
	#     bvr:   Breitung's variance ratio test.
	#     hurst: Hurst exponent as calculated by aggvarFit{fArma}.
	#
	# Returns an object of type htest representing the results of the hypothesis test.
	# In all cases, a low p.value is interpreted as evidence that the null hypothesis
	# of a unit root should be rejected.  In other words, a low p.value is evidence
    # for cointegration.
    #
    # For all of the tests except adfraw, the p-values have been re-calibrated.

	test.method <- match.arg(test.method)
    DNAME <- deparse(substitute(R))
	METHOD <- sprintf("Unit root test (%s) of residuals in Engle Granger procedure", test.method)

#    if (test.method %in% c("jo-e", "jo-t", "ers-p", "ers-d", "sp-r")) require(urca)
#    if (test.method == "hurst") require(fArma)

	if (test.method == "adfraw") {
		adf <- suppressWarnings(adf.test(R, "stationary"))
		STAT <- adf$statistic
		PVAL <- adf$p.value
		URTEST <- adf$method
	} else if (test.method == "adf") {
		adf <- suppressWarnings(adf.test(R, "stationary"))
		STAT <- adf$statistic
		PVAL <- quantile_table_interpolate(egc_adf_qtab, length(R), STAT)
		URTEST <- adf$method
	} else if (test.method == "pgff") {
		STAT <- pgff_rho_ws (R, detrend=TRUE)
		names(STAT) <- "rho_ws"
		PVAL <- quantile_table_interpolate(egc_pgff_qtab, length(R), STAT)
	    URTEST <- "Pantula, Gonzales-Farias and Fuller Unit Root Test"	
	} else if (test.method == "pp") {
		pp <- suppressWarnings(pp.test(R))
		STAT <- pp$statistic
		PVAL <- quantile_table_interpolate(egc_pp_qtab, length(R), STAT)
		URTEST <- pp$method	
	} else if (test.method == "ers-p") {
		ers <- ur.ers(R, type="P-test", model="constant")
		STAT <- ers@teststat
		names(STAT) <- "ERS P-test"
		PVAL <- quantile_table_interpolate(egc_ersp_qtab, length(R), STAT)
		URTEST <- ers@test.name		
	} else if (test.method == "ers-d") {
		ers <- ur.ers(R, type="DF-GLS", model="constant")
		STAT <- ers@teststat
		names(STAT) <- "ERS DF-GLS"
		PVAL <- quantile_table_interpolate(egc_ersd_qtab, length(R), STAT)
		URTEST <- ers@test.name		
	} else if (test.method == "sp-r") {
		sp <- ur.sp(R, type="rho")
		STAT <- sp@teststat
		names(STAT) <- "rho"
		PVAL <- quantile_table_interpolate(egc_spr_qtab, length(R), STAT)
		URTEST <- sp@test.name		
	} else if (test.method == "bvr") {
		STAT <- bvr_rho(R)
		names(STAT) <- "rho"
		PVAL <- quantile_table_interpolate(egc_bvr_qtab, length(R), STAT)
		URTEST <- "Breitung Variance Ratio Test"
	} else if (test.method == "hurst") {
		h <- aggvarFit(R)
		STAT <- h@hurst$H
		names(STAT) <- "H"
		PVAL <- quantile_table_interpolate(egc_hurst_qtab, length(R), STAT)
		URTEST <- h@title	
	} else {
		stop ("Unit root test.method ", test.method, " not implemented")
	}
		
    structure(list(statistic = STAT, alternative = "cointegrated", 
        p.value = PVAL, method = METHOD, urtest = URTEST, data.name = DNAME), 
        class = "htest")    
}

egc_test_power <- function (test.method=egcm.urtests(), 
	rho=0.95, n=250, nrep=10000, p.value=0.05) {
	# Uses simulation to estimate the power of a cointegration test. 
	# Power is defined to be the probability that the null hypothesis is rejected
	# when the null hypothesis is false.  Power is parameterized in terms of the
	# mean reversion coefficient rho, the sample size n, and the p-value.
	# The power is computed by generating random cointegrated 
	# pairs, and then counting the number of such pairs that are 
	# identified as cointegrated.

	test.method <- match.arg(test.method)
	pvalues <- replicate(nrep, test.egcm(rcoint(n, rho=rho), test.method)$p.value)
	sum(pvalues <= p.value) / length(pvalues)
}

egc_print_power_table <- function (pt) {
    cat(sprintf("%6s", ""))
    sapply(colnames(pt), function(cn) cat(sprintf(" %7s", cn)))
    cat("\n")
    for (i in 1:nrow(pt)) {
        cat(sprintf("%6s", rownames(pt)[i]))
        sapply(pt[i,], function(v) cat(sprintf(" %7.2f", v)))
        cat("\n")
    }
}

egc_test_power_table <- function (test.method=egcm.urtests(), 
	nrep=4000, p.value=0.05,
#	rho=c(0.995, 0.99, 0.98, 0.97, 0.96, 0.95, 0.90),
    rho=c(0.80, 0.90, 0.92, 0.94, 0.96, 0.98),
	n=c(60, 125, 250, 500, 1000)) {
	# Constructs a table of power estimates for realistic values of rho
	test.method <- match.arg(test.method)
	do.row <- function(nv) sapply(rho, function(r) egc_test_power(test.method, r, nv, nrep, p.value))
#    require(parallel)
	pt <- do.call("rbind", mclapply(n, function(nv) do.row(nv)))
	pt <- as.data.frame(pt)
	rownames(pt) <- n
	colnames(pt) <- rho
    egc_print_power_table(pt)
	invisible(pt)
}

egcm_power_comparison_table <- function (tests=egcm.urtests(), nrep=10000, p.value=0.05,
	rho=seq(0.85,0.99, by=0.01), n=250) {
	# Creates a table comparing the powers of various unit root tests for a fixed
	# sample size over a range of values of rho.

	do_test_method <- function (tm) {
		powers <- sapply(rho, function(r) egc_test_power(test.method=tm, rho=r, 
			n=n, nrep=nrep, p.value=p.value))
		df <- data.frame(Test=tm, Rho=rho, Power=powers)
		rownames(df) <- NULL
		df
	}
    
#	require(parallel)
	power_tab <- do.call("rbind", mclapply(tests, do_test_method))
	
	attr(power_tab, "n") <- n
	attr(power_tab, "nrep") <- nrep
	attr(power_tab, "p.value") <- p.value
	
	power_tab
}
	
egcm_power_graph <- function (power_tab = egcm_power_comparison_table()) {
	# Graphs the power estimates that were previously created with a call to
	# egcm_power_comparison_table()

    Rho <- NULL
    Power <- NULL
    Test <- NULL
    
	ggplot(power_tab, aes(x=Rho, y=Power, colour=Test)) + geom_line() +
		ggtitle(sprintf("Comparison of Powers of Cointegration Tests\nSample Size = %d", attr(power_tab, "n"))) 
		
}

egc_test_specificity <- function (test.method=egcm.urtests(), 
	nrep=1000, seqlen=250, p.value=c(0.01,0.02,0.05,0.1,0.2,0.5)) {
	# Uses simulation to estimate the specificity of a cointegration test.
	# Specificity is defined to be the probability that the null hypothesis
	# is accepted when the null hypothesis is true.  If the test is
	# calibrated properly, this will be identical to (one minus) the p-value 
	# reported by the test.  Thus, the purpose of this function is to check that
	# critical values are being calculated accurately.
	#
	# Input values:
	#   test.method:  Type of unit root test to be performed.  See egc.test for details.
	#   nrep:         Number of repetitions to perform
	#   seqlen:       Length of each randomly generated cointegrated pair
	#   p.value:      List of critical values that are to be checked
	#
	# For each p.value in the list, returns the percentage of trials where
	# the reported p-value was less than or equal to the specified value.
	# In addition, returns the standard errors and the significance levels.
	# The significance level is the probability that an observation this
	# extreme or greater would be observed (e.g., the two-tailed significance).
	
	test.method <- match.arg(test.method)
	replicates <- replicate(nrep, test.egcm(rcoint(seqlen, rho=1), test.method=test.method)$p.value)
	observed <- sapply(p.value, function(p) sum(replicates <= p))
	freq <- observed / nrep
	expected <- p.value * nrep
	se <- sqrt(p.value * (1 - p.value)/nrep)
	deviation <- abs(observed - expected)
	significance <- sapply(1:length(p.value),
		function(i) (1 - pbinom(expected[i]+deviation[i], nrep, p.value[i])) +
					pbinom(expected[i]-deviation[i], nrep, p.value[i]))
	rmat <- rbind(freq, se, significance)
	colnames(rmat) <- p.value
	rownames(rmat) <- c("Observed", "Standard Error", "Significance")
	rmat
}

egc_quantiles <- function(test.method=egcm.urtests(), 
	sample_size=100, nrep=40000, 
	q=c(0.0001, 0.001, 0.01, 0.025, 0.05, 0.10, 0.20, 0.50, 0.80, 0.90, 0.95, 0.975, 0.99, 0.999, 0.9999)) {
	# Calculates quantiles of the unit root test statistic under the assumption rho=1.
	test.method <- match.arg(test.method)
#	qvals <- replicate(nrep, test.egcm(rcoint(sample_size, alpha=0, beta=0, rho=1), test.method)$statistic)
	qvals <- replicate(nrep, egcm(rcoint(sample_size, rho=1), urtest=test.method)$r.stat)
	quantile(qvals, q)
}

egc_quantile_table <- function(test.method=egcm.urtests(), 
	nrep=40000,
	q=c(0.0001, 0.01, 0.025, 0.05, 0.10, 0.20, 0.50, 0.80, 0.90, 0.95, 0.975, 0.99, 0.999, 0.9999),
	n=c(25, 50, 100, 250, 500, 750, 1000, 1250, 2500)) {
	# Calculates a table of quantile values by sample size of the egc.test function
	# under the assumption rho=1.
	test.method <- match.arg(test.method)
#    require(parallel)
	df <- do.call("cbind", mclapply(n, function(nv) c(nv, egc_quantiles(test.method, nv, nrep, q))))
	df <- as.data.frame(df)
	colnames(df) <- n
	df <- cbind(data.frame(quantile=c(NA,q)), df) 
	df
}

################################################
##### Engle Granger Cointegration Model    #####
################################################

egcm <- function (X, Y, na.action, log=FALSE, normalize=FALSE, debias=TRUE, robust=FALSE, include.const=TRUE,
	i1test=egcm.default.i1test(), urtest=egcm.default.urtest(), p.value=egcm.default.pvalue()
    ) {
	# Performs the two-step Engle Granger cointegration analysis on the price
	# series X and Y, and creates an object representing the results of the
	# analysis.
	#  
	# If X is the price series of the first security and Y is 
	# the price series of the second, then computes the fit:
	#
	#     Y = alpha + beta * X + R
	#     R_t = rho * R_{t-1} + eps_t
	#
	# If log is TRUE, then the price series are logged before the analysis is 
	# performed.  If Y is missing, then X is presumed to be a two-column
	# matrix, and X[,1] and X[,2] are used in place of X and Y.
	
	if (missing(Y)) {
		if (is.null(ncol(X)) || (ncol(X) != 2)) {
			stop("If Y is missing, X must be a two-column matrix or data.frame")
		}
		series_names <- colnames(X)
		S2 <- X[,2]
		S1 <- X[,1]
	} else if (is.character(Y)) {
        stop("Y must be a numeric vector.  Did you mean to use yegcm?")
    } else if (length(Y) != length(X)) {
        stop("X and Y must be numeric vectors of the same length")
    } else {
        S1 <- X
        S2 <- Y
		series_names <- c(colnames(S1), colnames(S2))
	}

    if (!is.logical(log)) stop("Parameter log should be of type logical but got ", log)
    if (!is.logical(normalize)) stop("Parameter normalize should be of type logical but got ", normalize)
    if (!is.logical(debias)) stop("Parameter debias should be of type logical but got ", debias)
    if (p.value < 0.001) stop("P-values less than 0.001 are not supported")
    if (!is.logical(robust)) stop("Parameter robust should be of type logical but got ", robust)
    
	i1test <- match.arg(i1test, egcm.i1tests())
	urtest <- match.arg(urtest, egcm.urtests())

    if (length(series_names) == 0) series_names <- c("X", "Y")
	if (is.null(series_names[1])) series_names[1] <- "X"
	if (is.null(series_names[2])) series_names[2] <- "Y"
    if (series_names[1] == series_names[2]) series_names <- c("X","Y")
	
	if (normalize) {
		S1 <- S1/as.numeric(coredata(S1[which(!is.na(S1))[1]]))
		S2 <- S2/as.numeric(coredata(S2[which(!is.na(S2))[1]]))
	}
	
	if (log) {
		S1 <- log(S1)
		S2 <- log(S2)
	}
	
    if (missing(na.action)) na.action <- options("na.action")
    if (is.list(na.action) && is.character(na.action[[1]])) na.action <- na.action[[1]]
    if (is.character(na.action)) na.action <- get(na.action)
    S12 <- cbind(S1, S2)
    S12.na_action <- na.action(S12)
    S1 <- S12.na_action[,1]
    S2 <- S12.na_action[,2]

	if (is(S1, "zoo")) {
		S1index <- index(S1)
	} else {
		S1index <- 1:length(S1)
	}
	
	S1 <- coredata(S1)
	S2 <- coredata(S2)

    if (robust && include.const) {
        L <- summary(rlm(S2~S1))
    	alpha <- coef(L)[1,1]
    	beta <- coef(L)[2,1]        
    } else if (robust) {
        L <- summary(rlm(S2~S1+0))
    	alpha <- 0
    	beta <- coef(L)[1,1]       
    } else if (include.const) {
        L <- summary(lm(S2~S1))
    	alpha <- coef(L)[1,1]
    	beta <- coef(L)[2,1]        
    } else {
        L <- summary(lm(S2~S1+0))
    	alpha <- 0
    	beta <- coef(L)[1,1]    
    }
    
    N <- length(L$residuals)
    R <- L$residuals
    FR <- R[2:N]
    BR <- L$residuals[1:(N-1)]
    if (!robust) {
    	LR <- summary(lm(FR ~ BR + 0))
    	rho.raw <- coef(LR)[1,1]
    } else {
    	LR <- summary(rlm(FR ~ BR + 0))
    	rho.raw <- coef(LR)[1]
    }        

	if (debias) {
		rho <- debias_rho(length(FR), rho.raw)
	} else {
		rho <- rho.raw
	}
    eps <- FR - rho * BR
	
    if (include.const) {
        #	alpha.se <- coef(L)[1,2]
    	# The following works well in simulated data, but I am not sure if 
    	# it is correct.  In any event, simply taking the standard error
    	# generated by lm() is clearly incorrect.
    	alpha.se <- sqrt(coef(L)[1,2]^2 + var(eps)/4)
	
    	# Engle and Granger show that the convergence of beta is super-consistent,
    	# however simulation shows that the standard errors computed by lm() seem
    	# to be too small.  Perhaps a correction based on the var(eps) needs to
    	# be included here as well?
    	beta.se <- coef(L)[2,2]
    } else {
        alpha.se <- 0
        beta.se <- coef(L)[1,2]
    }
	
	rho.se <- coef(LR)[1,2]
#	rho.se <- rho_se(length(eps), rho) * sd(R)

	i1testfunc <- function(X) {
		switch(i1test,
			adf = suppressWarnings(adf.test(X)),
			pgff = pgff.test(X),
			bvr = bvr.test(X),
			pp = pp.test(X))
	}
	
	S1i1.test <- i1testfunc(S1)
	S2i1.test <- i1testfunc(S2)

	if (urtest %in% c("jo-t", "jo-e")) {
		R.test <- test.egcm(cbind(S1, S2), urtest)
	} else {
		R.test <- egc.residuals.test(R, urtest)
	}
	
	lb <- Box.test(eps, type="Ljung-Box")
	
	structure(list(S1=S1, S2=S2, 
        residuals=R, 
        innovations=eps, 
		series_names=series_names,
		index=S1index,
		i1test=i1test,
		urtest=urtest,
        pvalue=p.value,
		log=log,
		alpha=alpha,
		alpha.se = alpha.se,
		beta = beta,
		beta.se = beta.se,
		rho = rho,
		rho.raw=rho.raw,
		rho.se = rho.se,
		s1.i1.stat=S1i1.test$statistic,
		s1.i1.p=S1i1.test$p.value,
		s2.i1.stat=S2i1.test$statistic,
		s2.i1.p=S2i1.test$p.value,
		r.stat=R.test$statistic,
		r.p=R.test$p.value,
		eps.ljungbox.stat = lb$statistic,
		eps.ljungbox.p = lb$p.value,
        s1.dsd = sd(diff(S1)),
        s2.dsd = sd(diff(S2)),
        residuals.sd = sd(R),
		eps.sd = sd(eps)),
		class = "egcm")
}


print.egcm <- function (x, ...) {
    E <- x
    sign_string <- function (x) ifelse(x < 0, "-", "+")
    cat(sprintf("%s[i] = %8.4f %s[i] %s %8.4f + R[i], R[i] = %8.4f R[i-1] + eps[i], eps ~ N(0,%8.4f^2)\n",
        E$series_names[2], E$beta, E$series_names[1], 
        sign_string(E$alpha), abs(E$alpha),
        E$rho, E$eps.sd))

    pad <- function (s) paste0(rep(" ", nchar(s)), collapse="")
    alpha.se.str <- sprintf("(%.4f)", E$alpha.se)
    beta.se.str <- sprintf("(%.4f)", E$beta.se)
    rho.se.str <- sprintf("(%.4f)", E$rho.se)
    cat(sprintf("%s %14s%s %14s %23s\n",
        pad(E$series_names[2]), beta.se.str, pad(E$series_names[1]), alpha.se.str, rho.se.str))

#    cat(sprintf("R[LAST] = %.4f\n", last(E$residuals))) 

    N <- length(E$residuals)
#    cat(sprintf("\n%s[N] = %8.4f %s[N] %s %8.4f %s %8.4f\n",
#        E$series_names[2], E$beta, E$series_names[1], 
#        sign_string(E$alpha), abs(E$alpha),
#        sign_string(E$residuals[N]), abs(E$residuals[N])))
    
    cat(sprintf("\nR[%s] = %.4f (t = %.3f)\n", 
        format(E$index[N]),
        E$residuals[N], 
        E$residuals[N]/sd(E$residuals)))
    
    warnings <- c()
    if (E$s1.i1.p < E$pvalue) {
        warnings <- c(warnings, sprintf("%s does not seem to be integrated.", E$series_names[1]))
    } 
    if (E$s2.i1.p < E$pvalue) {
        warnings <- c(warnings, sprintf("%s does not seem to be integrated.", E$series_names[2]))
    } 
    if (!is.cointegrated(E)) {
        warnings <- c(warnings, sprintf("%s and %s do not appear to be cointegrated.",
            E$series_names[1], E$series_names[2]))
    } else if (!is.ar1(E)) {
        warnings <- c(warnings, "The series seem cointegrated but the residuals are not AR(1).")
    }
    
    if (length(warnings) > 0) {
        cat(sprintf("\nWARNING: %s\n", paste0(warnings, collapse=" ")))
    }
    
}

plot.egcm <- function (x, ...) plot.egcm.internal(x, ...)

plot.egcm.internal <- function (E, series_names=NULL, ...) {
    # The following definitions are made to prevent 'R CMD check' from barfing
    Date <- NULL
    Value <- NULL
    Facet <- NULL
    msd <- NULL
    # end of debarfication
    
	grid.newpage()
	pushViewport(viewport(layout=grid.layout(10,1)))
	
	if (missing(series_names)) {
		if (is.null(E$series_names)) {
			series_names <- c("S1", "S2")
		} else {
			series_names <- E$series_names
		}
	}
	
	s1.df <- data.frame(Date=E$index, Value=as.numeric(E$S1), Facet=series_names[1], Series="Price")
	s2.df <- data.frame(Date=E$index, Value=as.numeric(E$S2), Facet=series_names[2], Series="Price")
	s1scaled <- as.numeric(E$alpha + E$beta * E$S1)
    sign_string <- function (x) ifelse(x < 0, "-", "+")
	s1scaled_desc <- sprintf("%.2f %s %.2f * %s", E$alpha, sign_string(E$beta), abs(E$beta), series_names[1])
	s1scaled.df <- data.frame(Date=E$index, Value=s1scaled, Facet=s1scaled_desc,
		Series="Price")
	n <- length(E$index)
	R.df <- data.frame(Date=E$index, Value=E$residuals, Facet="R", 
		Series="Residuals")
	eps.df <- data.frame(Date=E$index[2:n], Value=E$innovations, Facet="Eps", 
		Series="Innovations")
	
	df1 <- rbind(s1scaled.df, s2.df)
	if (E$log) {
		ylabel <- "Log Price"
	} else {
		ylabel <- "Price"
	}
	p1 <- ggplot(df1, aes(x=Date, y=Value, colour=Facet)) + geom_line() +
		ylab(ylabel) + xlab("") +
		theme(legend.position="top") +
		scale_colour_discrete(name="") +
		ggtitle ("Price Series")
        if(!is.cointegrated(E)) {
            x <- min(df1$Date)
            ymin <- min(df1$Value)
            ymax <- max(df1$Value)
            y <- ymin + 0.96 * (ymax - ymin)
            p1 <- p1 + annotate("text", x=x, y=y, label="Not cointegrated", colour="red", hjust=0, size=3)
        }
	print(p1, vp=viewport(layout.pos.row=1:4, layout.pos.col=1))

	sdR <- sd(E$residuals)
	hlines = data.frame(Value=c(2 * sdR, sdR, -sdR, -2 * sdR),
		Facet=c("two", "one", "one", "two"))
	p2 <- ggplot(R.df, aes(x=Date, y=Value)) + geom_line() +
		ggtitle ("Residual Series") + ylab("Differential") + xlab("") +
		geom_hline(data=hlines, aes(yintercept=Value, colour=Facet), linetype="dashed") +
        guides(colour=FALSE)

    sdstr <- sprintf("sd(R) = %.2f", E$residuals.sd)
    x <- min(R.df$Date)
    ymin <- min(R.df$Value)
    ymax <- max(R.df$Value)
    y <- ymin + 0.96 * (ymax - ymin)
    p2 <- p2 + annotate("text", x=x, y=y, label=sdstr, hjust=0, size=3)

	print(p2, vp=viewport(layout.pos.row=5:7, layout.pos.col=1))
	
	eps.df$sd <- rollapply(eps.df$Value, 20, sd, align="right", fill=NA)
	eps.df$msd <- -eps.df$sd
	p3 <- ggplot(eps.df, aes(x=Date, y=Value)) + geom_line() +
		ggtitle ("Innovations") + ylab("Epsilon") + xlab("") +
		geom_ribbon(aes(ymin=msd, ymax=sd), alpha=0.3, na.rm=TRUE)
	print(p3, vp=viewport(layout.pos.row=8:10, layout.pos.col=1))
}

as.data.frame.egcm <- function (x, row.names, optional, ...) {
    # The purpose of the following is to prevent 'R CMD check' from
    # barfing over the fact that the method signature for as.data.frame
    # does not match the signature of as.data.frame.egcm.internal.
    m <- match.call()
    args <- list(...)
    args$E <- x
    args$row.names <- m$row.names
    do.call(as.data.frame.egcm.internal, args, envir=sys.frame(sys.parent(2)))
}

as.data.frame.egcm.internal <- function (E, row.names=NULL, i1tests=c(), urtests=c()) {
	df <- data.frame(series1=E$series_names[1],
		series2=E$series_names[2],
		log=E$log,
		i1test=E$i1test,
		urtest=E$urtest,
		alpha=E$alpha, 
		alpha.se=E$alpha.se, 
	    beta=E$beta, 
	    beta.se=E$beta.se,
	    rho=E$rho, 
		rho.raw=E$rho.raw,
	    rho.se=E$rho.se, 
		s1.i1.stat=E$s1.i1.stat,
		s1.i1.p=E$s1.i1.p,
		s2.i1.stat=E$s2.i1.stat,
		s2.i1.p=E$s2.i1.p,
		r.stat=E$r.stat,
		r.p=E$r.p,
		eps.ljungbox.stat=E$eps.ljungbox.stat,
		eps.ljungbox.p=E$eps.ljungbox.p,
        s1.dsd=E$s1.dsd,
        s2.dsd=E$s2.dsd,
        residuals.sd = E$residuals.sd,	
	    eps.sd=E$eps.sd)
	rownames(df) <- row.names
	
	i1testfunc <- function(X, testtype) {
		switch(testtype,
			adf = suppressWarnings(adf.test(X)),
			pgff = pgff.test(X),
			bvr = bvr.test(X),
			pp = pp.test(X))
	}

	for (i1t in i1tests) {
		i1r <- i1testfunc(E$S1, i1t)
		i2r <- i1testfunc(E$S2, i1t)
		dfcols <- paste0(c("s1.", "s1.", "s2.", "s2."), i1t, c(".stat", ".p", ".stat", ".p"))
		df[,dfcols] <- c(i1r$statistic, i1r$p.value, i2r$statistic, i2r$p.value)
	}
	
	for (ut in urtests) {
		ur <- test.egcm(E, ut)
		dfcols <- paste0("r.", ut, c(".stat", ".p"))
		df[,dfcols] <- c(ur$statistic, ur$p.value)
	}
	
	df
}

summary.egcm <- function (object, ...) {
	structure(list(EGCM=object), class="summary.egcm")	
}

print.summary.egcm <- function (x, ...) {
    ES <- x
	print(ES$EGCM)
	cat("\nUnit Root Tests of Residuals\n")
	cat(sprintf("%-50s %10s %10s\n", "", "Statistic", "p-value"))
	do.test <- function(testname, testdesc) {
		result <- test.egcm(ES$EGCM, testname)
		cat(sprintf("  %-48s %10.3f %10.5f\n", testdesc, result$statistic, result$p.value))
	}
	do.test("adf",   "Augmented Dickey Fuller (ADF)")
	do.test("pp",    "Phillips-Perron (PP)")
	do.test("pgff",  "Pantula, Gonzales-Farias and Fuller (PGFF)")
	if (exists("ur.ers")) do.test("ers-d", "Elliott, Rothenberg and Stock DF-GLS (ERSD)")
	if (exists("ca.jo"))  do.test("jo-t",  "Johansen's Trace Test (JOT)")
	if (exists("ur.sp"))  do.test("sp-r",  "Schmidt and Phillips Rho (SPR)")
	
	cat("\n")
	cat("Variances\n")
	vs1 <- sprintf("SD(diff(%s))", ES$EGCM$series_names[1])
	vs2 <- sprintf("SD(diff(%s))", ES$EGCM$series_names[2])
	cat(sprintf("  %-20s = %10.6f\n", vs1, sd(diff(ES$EGCM$S1))))
	cat(sprintf("  %-20s = %10.6f\n", vs2, sd(diff(ES$EGCM$S2))))
	cat(sprintf("  %-20s = %10.6f\n", "SD(diff(residuals))", sd(diff(ES$EGCM$residuals))))
	cat(sprintf("  %-20s = %10.6f\n", "SD(residuals)", sd(ES$EGCM$residuals)))
	cat(sprintf("  %-20s = %10.6f\n", "SD(innovations)", sd(ES$EGCM$innovations)))
	
	if (ES$EGCM$rho < 1) {
		half_life <- log(0.5) / log(ES$EGCM$rho)
		cat(sprintf("\nHalf life       = %10.6f\n", half_life))
	} else {
		half_life <- Inf
		cat(sprintf("\nHalf life       = Infinite\n"))
	}
	Rlast = ES$EGCM$residuals[length(ES$EGCM$residuals)]
	Rlast_t = Rlast / sd(ES$EGCM$residuals)
	cat(sprintf("R[last]         = %10.6f (t=%.2f)\n", Rlast, Rlast_t))
	
	if (!ES$EGCM$log && include.investment.scenarios()) {
		cat("\nInvestment scenarios\n")
		isc <- data.frame(open_sd=c(2.0, 2.0, 1.5, 1.0), close_sd=c(1.0, 0.5, 0.5, 0.5))
		sdr <- sd(ES$EGCM$residuals)
		n <- length(ES$EGCM$S1)
		capital_req <- ES$EGCM$S1[n]*ES$EGCM$beta + ES$EGCM$S2[n]
		if (ES$EGCM$alpha < 0) {
			capital_req <- capital_req + ES$EGCM$alpha
		}
		cat(sprintf("  Capital requirement = %8.2f\n\n", capital_req))
		cat(sprintf(" %10s %6s  %10s %6s %10s %20s %15s %12s\n",
			"Open", "(SD)", "Close", "(SD)", "E(Gain)", "E(Return on Capital)", 
			"E(Holding Time)", "Ann Return"))
		for (i in 1:nrow(isc)) {
			open_r <- isc$open_sd[i] * sdr
			open_sd_str <- sprintf("(%.1f)", isc$open_sd[i])
			close_r <- isc$close_sd[i] * sdr
			close_sd_str <- sprintf("(%.1f)", isc$close_sd[i])
			egain <- open_r - close_r
			eroc <- egain / capital_req
			eht <- half_life * log(open_r/close_r)/log(2)
			ann_return <- 100.0 * (exp(log(1+eroc) * 252/eht) - 1)
			cat(sprintf(" %10.3f %6s  %10.3f %6s %10.3f %20.5f %15.3f %12.2f%%\n",
				open_r, open_sd_str, close_r, close_sd_str,
				egain, eroc, eht, ann_return))
		}
	}
}

sim.egcm <- function (E, nsteps, X0, Y0) {
    # Given a cointegration model, simulates it for nsteps,
    # starting with initial values (X0, Y0).  If X0 or Y0
    # is missing, then simulates from the last values used
    # in generating the model E.  Returns a two-column
    # data.frame, where the first column contains the simulated
    # values of X0, and the second column contains the simulated
    # values of Y0.  
    if (missing(X0) || missing(Y0)) {
        X0 <- E$S1[length(E$S1)]
        Y0 <- E$S2[length(E$S2)]
    }
    s <- rcoint (nsteps, E$alpha, E$beta, E$rho, E$s1.dsd, E$eps.sd, X0, Y0)
    colnames(s) <- E$series_name
    rownames(s) <- NULL
    s
}

is.cointegrated <- function (E) {
	# Returns TRUE if the egcm model E appears to be cointegrated.
	
	if (!is(E,"egcm")) {
		stop("Parameter E must be of type egcm")
	}

	S1.I1 <- E$s1.i1.p > E$pvalue
	S2.I1 <- E$s2.i1.p > E$pvalue
	R.I0 <- E$r.p < E$pvalue
	
	S1.I1 && S2.I1 && R.I0
}

is.ar1 <- function(E) {
	# Returns TRUE if the residuals in the egcm model E are adequately
	# fit by an AR(1) model.
	
	if (!is(E,"egcm")) {
		stop("Parameter E must be of type egcm")
	}
	
	return (E$eps.ljungbox.p > E$pvalue)
}

residuals.egcm <- function(object, ...) {
    # Returns the residual series R[i] associated to an egcm object
    
    E <- object
	if (!is(E,"egcm")) {
		stop("Parameter E must be of type egcm")
	}
    return (E$residuals)
}

################################################
##### De-biasing of Rho                    #####
################################################

rho_bias_table <- function(N=c(60,125,250,500,750,1000,2500), 
	RHO=c(0.9, 0.95, 0.96, 0.97, 0.98, 0.99, 1),
    nsamples=1000,
	print=TRUE,
	report.sd=FALSE) {
	# Calculates the bias of an estimator RHO_EST of rho in the Engle-Granger model
	#   Y[t] = alpha + beta * X[t] + R[t]
	#   R[t] = rho * R[t-1] + eps
	# where eps is NID(0,1).  
	#
	# N specifies the various vector lengths that are to be considered,
	# and RHO specifies the actual values of RHO that are to be used in
	# the simulation.

	rho_hat_sample <- function (sample_size, rho) {
		X <- rcoint(sample_size, rho=rho)
		S2 <- X[,2]
		S1 <- X[,1]
		L <- summary(lm(S2 ~ S1))
		R <- L$residuals
		BR <- c(R[2:length(R)], NA)
		LR <- summary(lm(R ~ BR + 0))
		coef(LR)[1,1]
	}

	
	bias <- function(n, r) {
#		r_hat <- replicate(nsamples, rho_hat_sample(n, r))
		pvec_rho_hat_sample <- function (v) {
			replicate(length(v), rho_hat_sample(n, r))
		}
		r_hat <- parallel::pvec(1:nsamples, pvec_rho_hat_sample)
		c(mean(r_hat), sd(r_hat))
	}

	M <- matrix(0, ncol=length(N), nrow=length(RHO))
	
	if (print) {
		cat(sprintf("%4s ", "N"))
		for (r in RHO) {
			cat(sprintf("%7.4f %9s ", r, ""))
		}
		cat("\n")
	}
	for (i in index(N)) {
		n <- N[i]
		if (print) cat(sprintf("%4d ", n))
		for (j in index(RHO)) {
			r <- RHO[j]
			b <- bias(n,r)
			if (report.sd) {
				M[j,i] <- b[2] * sqrt(n)
			} else {
				M[j,i] <- b[1]				
			}
			if (print) cat(sprintf("%7.4f (%7.4f) ", b[1], b[2]))
		}
		if (print) cat("\n")
	}
	
	M <- rbind(N, M)
	M <- cbind(c(NA, RHO), M)
	
	colnames(M) <- c("quantile", N)
	rownames(M) <- c("", RHO)
	
	M
}

qtab_to_ltab <- function (tab=rho_bias_qtab) {
	do.call("rbind", lapply(2:ncol(tab), 
		function(i) {
			slm <- summary(lm(tab[2:nrow(tab),1] ~ tab[2:nrow(tab), i]))
			n = tab[1,i]
			a0 = coef(slm)[1,1]
			a1 = coef(slm)[2,1]
			data.frame(n, a0, a1)
		}))
}

generate_rho_bias_table <- function() {
	qtab <- rho_bias_table(N=c(25, 50, 100, 200, 400, 800, 1200, 1600, 2000), 
	  RHO=c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9),
      nsamples=10000, print=FALSE)
	qtab_to_ltab(qtab)
}

debias_rho <- function (sample_size, rho) {
	# Given an estimate rho of the mean reversion parameter rho obtained through
	# the Engle Granger procedure, calculates a de-biased estimate of rho.
	#
	# Several different approaches have been tried for debiasing rho.  The
	# first approach was to compute a table of values indexed by sample_size
	# and rho.  Each entry in the table contains the mean of the estimated
	# value rho_hat for that pair (sample_size, rho).  The table is then
	# inverted using linear interpolation to obtain a debiased estimate of
	# rho_hat.  This approach was found to yield an improvement, but the
	# estimates obtained in this way were still biased.  
	#
	# It was then observed that for a given sample_size, the entries in
	# the debias table are nearly perfectly linear in rho.  Based upon this,
	# the debias table was condensed down to a set of linear relations,
	# indexed by sample_size.  Interpolation is used for sample_sizes that
	# are intermediate with respect to those represented in the table.
	#
	# Upon inspection of the debiasing table, it was found that the following
	# relation seems to offer a reasonably good de-biasing estimate:
	#    rho_hat_debiased = rho_hat * (1 + 6/sample_size)
	#
	# Another approach that was tried was to choose random cointegrated
	# vectors and then to use the Engle Granger procedure to estimate 
	# rho_hat.  A table of values (rho, rho_hat) was constructed in this
	# way, and then coefficients c_0 and c_1 were obtained which gave
	# the best fit to rho_hat = c_0 + c_1 rho.  These values were then
	# inverted to obtain the debiasing relation, e.g.,
	#    rho_hat_debiased = -(c_0/c_1) + (1 / c_1) rho_hat
	# This last approach seems to give the best fit, and so it is the
	# method that is used.
	
#	return(rho * (1 + 6/sample_size))
#	return(quantile_table_interpolate(rho_bias_qtab, sample_size, rho))
	
	if (rho < 0) {
		-debias_rho (sample_size, -rho)
	} 

	# For a given sample size, the bias is well described as a linear
	# function of rho.
	i <- findInterval(sample_size, rho_bias_ltab$n)
	if (i == 0) {
		y <- rho_bias_ltab$c0[1] + rho_bias_ltab$c1[1] * rho
	} else if (i < nrow(rho_bias_ltab)) {
		y1 <- rho_bias_ltab$c0[i] + rho_bias_ltab$c1[i] * rho
		y2 <- rho_bias_ltab$c0[i+1] + rho_bias_ltab$c1[i+1] * rho
		n1 <- rho_bias_ltab$n[i]
		n2 <- rho_bias_ltab$n[i+1]
		w1 <- (n2 - sample_size) / (n2 - n1)
		w2 <- 1 - w1
		y <- y1 * w1 + y2 * w2
		y <- min(y, 1)
	} else {
		y1 <- rho_bias_ltab$c0[i] + rho_bias_ltab$c1[i] * rho
		w1 <- exp(1-sample_size/rho_bias_ltab$n[i])
		w2 <- 1 - w1
		y <- y1 * w1 + rho * w2
		y <- min(y, 1)
	}
	
	y
}

coint_bias_table <- function (sample_size=c(60, 125, 250, 500, 750, 1000), 
  nsamples=10000, debias=TRUE, intercept=TRUE, debug=FALSE, reverse=FALSE) {
	# Generates nsamples random pairs of cointegrated vectors of length sample_size, 
	# and uses egcm to estimate alpha, beta and rho.  Then, calculates the regressions:
	#    alpha_hat = a_0 + a_1 * alpha
	#    beta_hat = b_0 + b_1 * beta
	#    rho_hat = c_0 + c_1 * rho
	# If the estimates are unbiased, then we should find that
	#    a_0 = 0, a_1 = 1,
	#    b_0 = 0, b_1 = 1, and
	#    c_0 = 0, c_1 = 1.
	# The output can be inspected to see whether or not these relations hold.
	#
	# Parameters are as follows:
	#    sample_size:    List of sample sizes for which biases are to be computed.
	#    nsamples:       Number of iterations to perform for each sample size.
	#    debias:         If TRUE, then egcm attempts to compute debiased estimates
	#                    of rho, and the biases of these debiased estimates are reported.
	#    intercept:      If TRUE, then a constant term is included in the bias estimates
	#    reverse:        If TRUE, then the linear models are inverted, e.g., the following
	#                    models are computed:
	#                         alpha = a_0 + a_1 * alpha_hat, etc.
	#    debug:          If TRUE, then the debugger is invoked after the random samples
	#                    are computed for each sample_size.
	#
	# Returns a data.frame with the bias estimates.
	
	do_rep <- function (ss) {
		alpha <- runif(1,-10,10) 
		beta <- runif(1,-10,10)
		rho <- runif(1,0,1)	
		X <- rcoint(ss, alpha, beta, rho)
		e <- egcm(X, debias=debias)
		data.frame(alpha_hat=e$alpha, alpha=alpha, 
			       beta_hat=e$beta, beta=beta,
			       rho_hat=e$rho, rho=rho, rho.raw=e$rho.raw)
	}
	
	recenter <- function (cvec, loc=1, ndf=30) {
		rvec <- cvec
		rvec[3] <- (rvec[1]-loc)/rvec[2]
		rvec[4] <- pt(-abs(rvec[3]), ndf)
		rvec
	}
	
	sample_summary <- function (data, ssize) {
		coef.names <- c("a0", "a1", "b0", "b1", "c0", "c1", "c0raw", "c1raw")

		if (reverse) {
			s.alpha <- summary(lm(alpha ~ alpha_hat, data))
			s.beta <- summary(lm(beta ~ beta_hat, data))
			s.rho <- summary(lm(rho ~ rho_hat, data))
			s.rho.raw <- summary(lm(rho ~ rho.raw, data))
		
			sm <- rbind(coef(s.alpha)[1,],
				        recenter(coef(s.alpha)[2,], ndf=s.alpha$df[2]),
						coef(s.beta)[1,],
						recenter(coef(s.beta)[2,], ndf=s.beta$df[2]),
						coef(s.rho)[1,],
						recenter(coef(s.rho)[2,], ndf=s.rho$df[2]),
						coef(s.rho.raw)[1,],
						recenter(coef(s.rho.raw)[2,], ndf=s.rho.raw$df[2]))
						
		} else if (intercept) {
			s.alpha <- summary(lm(alpha_hat ~ alpha, data))
			s.beta <- summary(lm(beta_hat ~ beta, data))
			s.rho <- summary(lm(rho_hat ~ rho, data))
			s.rho.raw <- summary(lm(rho.raw ~ rho, data))
		
			sm <- rbind(coef(s.alpha)[1,],
				        recenter(coef(s.alpha)[2,], ndf=s.alpha$df[2]),
						coef(s.beta)[1,],
						recenter(coef(s.beta)[2,], ndf=s.beta$df[2]),
						coef(s.rho)[1,],
						recenter(coef(s.rho)[2,], ndf=s.rho$df[2]),
						coef(s.rho.raw)[1,],
						recenter(coef(s.rho.raw)[2,], ndf=s.rho.raw$df[2]))
						
		} else {
			s.alpha <- summary(lm(alpha ~ alpha_hat + 0, data))
			s.beta <- summary(lm(beta ~ beta_hat + 0, data))
			s.rho <- summary(lm(rho ~ rho_hat + 0, data))
			s.rho.raw <- summary(lm(rho.raw ~ rho + 0, data))
		
			sm <- rbind(0,
				        recenter(coef(s.alpha)[1,], ndf=s.alpha$df[2]),
						0,
						recenter(coef(s.beta)[1,], ndf=s.beta$df[2]),
						0,
						recenter(coef(s.rho)[1,], ndf=s.rho$df[2]),			
						0,
						recenter(coef(s.rho.raw)[2,], ndf=s.rho.raw$df[2]))
		}
				
		sm <- as.data.frame(sm)
		sm$coef <- coef.names
		sm$sample_size <- ssize
		sm <- sm[,c(6,5,1,2,3,4)]
		sm
	}
	if (debug) debug(sample_summary)
#	require(parallel)
	ssfunc <- function(ss) {
		data <- do.call("rbind", mclapply(1:nsamples, function(i) do_rep(ss)))
		sample_summary(data, ss)			
	}
	
	do.call("rbind", lapply(sample_size, ssfunc))
}

cbtab_to_rhobtab <- function (cbt) {
        n <- cbt$sample_size[cbt$coef == "c0raw"]
        c0r <- cbt$Estimate[cbt$coef == "c0raw"]
        c1r <- cbt$Estimate[cbt$coef == "c1raw"]
        c0 <- -c0r / c1r
        c1 <- 1.0 / c1r
        data.frame(n, c0, c1)
}
