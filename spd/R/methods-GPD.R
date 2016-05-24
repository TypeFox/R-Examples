#################################################################################
##
##   R package spd by Alexios Ghalanos Copyright (C) 2008-2013
##   This file is part of the R package spd.
##
##   The R package spd is free software: you can redistribute it and/or modify
##   it under the terms of the GNU General Public License as published by
##   the Free Software Foundation, either version 3 of the License, or
##   (at your option) any later version.
##
##   The R package spd is distributed in the hope that it will be useful,
##   but WITHOUT ANY WARRANTY; without even the implied warranty of
##   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##   GNU General Public License for more details.
##
##   Code Included and Modified from fExtremes package by Diethelm Wuertz
##   and other authors from which code in that package that was derived.
#################################################################################
# Developer Notes: included the gpd distribution from the fExtremes package so as not
# to require its dependencies (fBasics..etc)
gpdfit<-function(x, u = quantile(x, 0.95), type = c("mle", "pwm"),
		information = c("observed", "expected"), title = NULL, description = NULL, ...)
{
	UseMethod("gpdfit")
}

.gpdfit = function(x, u = quantile(x, 0.95), type = c("mle", "pwm"), 
		information = c("observed", "expected"), title = NULL, description = NULL, ...)
{   # A function implemented by Diethelm Wuertz and included/modified by Alexios
	# Ghalanos
	# Description:
	#   Fits a generalized Pareto model to excesses
	# Details:
	#   Returns an object of class "fGPDFIT" representing the fit of a
	#   generalized Pareto model to excesses over a high threshold
	# Notes:
	#   This is a wrapper to EVIR's 'gpd' function.
	call = match.call()
	type = match.arg(type)
	information = match.arg(information)
	# Check Type and Convert:
	X = x
	x = as.vector(x)
	N = length(x)
	# Compute Exceedances:
	exceedances = x[x > u]
	Names = as.character((1:N)[x > u])
	exceedances = as.vector(exceedances)
	names(exceedances) = Names
	# Estimate Parameters:
	if (type == "mle") {
		fit = .gpdmleFit(x, u, information)
		fit$llh = fit$fit$value
		fit$convergence = fit$fit$convergence
	} else if (type == "pwm") {
		fit = .gpdpwmFit(x, u)
		fit$llh = NA
		fit$convergence = NA
	}
	fit$p.less.thresh = fit$prob = 1 - length(x[x > u]) / length(x)
	fit$threshold = u
	fit$data = x
	# Compute Residuals:
	xi = fit$par.ests["xi"]
	beta = fit$par.ests["beta"]
	residuals = log(1 + (xi * (as.vector(exceedances)-u))/beta)/xi
	# Add title and description:
	if (is.null(title)) title = "GPD Parameter Estimation"
	if (is.null(description)) description = .description()
	# Compose Parameter List:
	parameter = list(u = u, type = type)
	if (information == "mle") parameter$information = information
	# Return Value:
	new("GPDFIT",
			call = call,
			method = c("gpd", type),
			parameter = parameter,
			data = list(x = X, exceedances = exceedances),
			fit = fit,
			residuals = residuals,
			title = title,
			description = description)
}

setMethod(f="gpdfit",definition=.gpdfit)

.gpdmleFit <-
    function (x, u, information = c("observed", "expected"), ...)
{
    # A Copy from Evir

    data = x
    threshold = u

    exceedances <- data[data > threshold]
    excess <- exceedances - threshold
    Nu <- length(excess)
    xbar <- mean(excess)

    s2 <- var(excess)
    xi0 <- -0.5 * (((xbar * xbar)/s2) - 1)
    beta0 <- 0.5 * xbar * (((xbar * xbar)/s2) + 1)
    theta <- c(xi0, beta0)

    negloglik <- function(theta, tmp)
    {
        xi <- theta[1]
        beta <- theta[2]
        cond1 <- beta <= 0
        cond2 <- (xi <= 0) && (max(tmp) > (-beta/xi))
        if (cond1 || cond2) {
            f <- 1e+06
        } else {
            y <- logb(1 + (xi * tmp)/beta)
            y <- y/xi
            f <- length(tmp) * logb(beta) + (1 + xi) * sum(y)
        }
        f
    }

    fit <- optim(theta, negloglik, hessian = TRUE, ..., tmp = excess)
    names(fit$par) = c("xi", "beta")

    if (fit$convergence) warning("\nOptimization may not have succeeded.")
    par.ests <- fit$par
    converged <- fit$convergence
    nllh.final <- fit$value

    information <- match.arg(information)
    if (information == "observed")
        varcov <- solve(fit$hessian)
    if (information == "expected") {
        one <- (1 + par.ests[1])^2/Nu
        two <- (2 * (1 + par.ests[1]) * par.ests[2]^2)/Nu
        cov <- -((1 + par.ests[1]) * par.ests[2])/Nu
        varcov <- matrix(c(one, cov, cov, two), 2)
    }

    par.ses <- sqrt(diag(varcov))
    names(par.ses) = c("xi", "beta")

    list(par.ests = par.ests, par.ses = par.ses, fit = fit, varcov = varcov)
}


#-------------------------------------------------------------------------------
.gpdpwmFit <- function (x, u)
{
    # A Copy from Evir

    data = x
    threshold = u

    data <- as.numeric(data)
    n <- length(data)
    exceedances <- data[data > threshold]
    excess <- exceedances - threshold
    Nu <- length(excess)
    xbar <- mean(excess)

    a0 <- xbar
    gamma <- -0.35
    delta <- 0
    pvec <- ((1:Nu) + gamma)/(Nu + delta)
    a1 <- mean(sort(excess) * (1 - pvec))
    xi <- 2 - a0/(a0 - 2 * a1)
    beta <- (2 * a0 * a1)/(a0 - 2 * a1)
    par.ests <- c(xi, beta)
    names(par.ests) = c("xi", "beta")

    denom <- Nu * (1 - 2 * xi) * (3 - 2 * xi)
    if (xi > 0.5) {
        denom <- NA
        warning("Asymptotic Standard Errors not available for PWM when xi>0.5.")
    }
    one <- (1 - xi) * (1 - xi + 2 * xi^2) * (2 - xi)^2
    two <- (7 - 18 * xi + 11 * xi^2 - 2 * xi^3) * beta^2
    cov <- beta * (2 - xi) * (2 - 6 * xi + 7 * xi^2 - 2 * xi^3)
    varcov <- matrix(c(one, cov, cov, two), 2)/denom
    information <- "expected"
    converged <- NA
    nllh.final <- NA
    par.ses <- sqrt(diag(varcov))
    names(par.ses) = c("xi", "beta")

    list(par.ests = par.ests, par.ses = par.ses, fit = NA, varcov = NA)
}


################################################################################
dgpd<-function(x, xi = 1, mu = 0, beta = 1, log = FALSE)
{   # A function written by Diethelm Wuertz
	# Description:
	#   Density for the Generalized Pareto DF
	# Transform:
	shape = xi
	location = mu
	scale = beta
	# Density:
	d = .dgpd(x, location, scale, shape, log)
	# Return Value:
	return(d)
}
#-------------------------------------------------------------------------------
pgpd<-function(q, xi = 1, mu = 0, beta = 1, lower.tail = TRUE)
{   # A function written by Diethelm Wuertz
	# Description:
	#   Probability for the Generalized Pareto DF
	# FUNCTION:
	# Transform:
	shape = xi
	location = mu
	scale = beta
	# Probability:
	p = .pgpd(q, location, scale, shape, lower.tail)
	# Return Value:
	return(p)
}

#-------------------------------------------------------------------------------
qgpd<-function(p, xi = 1, mu = 0, beta = 1, lower.tail = TRUE)
{   # A function written by Diethelm Wuertz
	# Description:
	#   Quantiles for the Generalized Pareto DF
	# FUNCTION:
	# Transform:
	shape = xi
	location = mu
	scale = beta
	# Quantiles:
	q = .qgpd(p, location, scale, shape, lower.tail)
	# Return Value:
	return(q)
}

#-------------------------------------------------------------------------------
rgpd<-function(n, xi = 1, mu = 0, beta = 1)
{   # A function written by Diethelm Wuertz
	# Description:
	#   Random variates for the Generalized Pareto DF
	# FUNCTION:
	# Transform:
	shape = xi
	location = mu
	scale = beta
	# Random Variates:
	r = .rgpd(n, location, scale, shape)
	# Return Value:
	return(r)
}

#-------------------------------------------------------------------------------
.dgpd<-function(x, location = 0, scale = 1, shape = 0, log = FALSE)
{
	stopifnot(min(scale) > 0)
	stopifnot(length(shape) == 1)
	# Density:
	d <- (x - location)/scale
	nn <- length(d)
	scale <- rep(scale, length.out = nn)
	index <- (d > 0 & ((1 + shape * d) > 0)) | is.na(d)
	if (shape == 0) {
		d[index] <- log(1/scale[index]) - d[index]
		d[!index] <- -Inf
	} else {
		d[index] <- log(1/scale[index]) - (1/shape+1)*log(1+shape*d[index])
		d[!index] <- -Inf
	}

	# Log:
	if (!log) d <- exp(d)
	# Return Value:
	return(d)
}

.pgpd<-function(q, location = 0, scale = 1, shape = 0, lower.tail = TRUE)
{
	stopifnot(min(scale) > 0)
	stopifnot(length(shape) == 1)
	# Probability:
	q <- pmax(q - location, 0)/scale
	if (shape == 0)
		p <- 1 - exp(-q)
	else {
		p <- pmax(1 + shape * q, 0)
		p <- 1 - p^(-1/shape)
	}
	# Lower Tail:
	if (!lower.tail) p <- 1 - p
	return(p)
}

.qgpd<-function(p, location = 0, scale = 1, shape = 0, lower.tail = TRUE)
{
	stopifnot(min(scale) > 0)
	stopifnot(length(shape) == 1)
	stopifnot(min(p, na.rm = TRUE) >= 0)
	stopifnot(max(p, na.rm = TRUE) <= 1)
	# Lower Tail:
	if (lower.tail)
		p <- 1 - p
	# Quantiles:
	if (shape == 0) {
		q = location - scale * log(p)
	} else {
		q = location + scale * (p^(-shape) - 1)/shape
	}
	# Return Value:
	return(q)
}

.rgpd<-function(n, location = 0, scale = 1, shape = 0)
{
	# could have just called qgpd
	stopifnot(min(scale) > 0)
	stopifnot(length(shape) == 1)
	# Random Variates:
	if (shape == 0) {
		r = location + scale * rexp(n)
	} else {
		r = location + scale * (runif(n)^(-shape) - 1)/shape
	}
	# Return Value:
	return(r)
}

###########################################################################
# Plot Functions for GDP

.plot.GPDFIT<-function(x, which = "ask", ...)
{   # A function implemented by Diethelm Wuertz, modified by Alexios Ghalanos
	# Description:
	#   Plot method for objects of class 'GPDFIT'
		.interactivegpdfitPlot(
			object = x,
			choices = c(
					"Excess Distribution",
					"Tail of Underlying Distribution",
					"Scatterplot of Residuals",
					"QQ-Plot of Residuals"),
			paste(".plot.gpdfit.", 1:4, sep = ""),
			which = which, ...)

	# Return Value:
	invisible(x)
}

.interactivegpdfitPlot = function(object, choices, plotFUN, which, ...)
{
# Some cecks:
	if (length(choices) != length(plotFUN))
		stop("Arguments choices and plotFUN must be of same length.")
	if (length(which) > length(choices))
		stop("Arguments which has incorrect length.")
	if (length(which) > length(plotFUN))
		stop("Arguments which has incorrect length.")
	# Plot:
	if (is.numeric(which)) {
		Which = rep(FALSE, times = length(choices))
		Which[which] = TRUE
	}

	if (which[1] == "all") {
		Which = rep(TRUE, times = length(choices))
	}

	if (which[1] == "ask") {
		.multgpdfitPlot(object, choices, plotFUN, ...)
	} else {
		for ( i in 1:length(choices) ) {
			FUN = match.fun(plotFUN[i])
			if (Which[i]) FUN(object)
		}
	}

	# Return Value:
	invisible(object)
}

.multgpdfitPlot = function (object, choices, ...)
{
	pick = 1
	while (pick > 0)
	{
		pick = menu(
				### choices = paste("plot:", choices)
				choices = paste("plot:",choices,sep=""),
				title = "\nMake a plot selection (or 0 to exit):")
		switch (pick,
				.plot.gpdfit.1(object,...),  .plot.gpdfit.2(object,...),  .plot.gpdfit.3(object,...),
				.plot.gpdfit.4(object,...))
	}
}
# ------------------------------------------------------------------------------
.plot.gpdfit.1<-function(x, labels = TRUE, ...)
{   # A function implemented by Diethelm Wuertz
	# and modified by Alexios Ghalanos
	# Description:
	#   Empirical Distribution Plot
	# Arguments:
	#   x - an object of class fGPDFIT
	#   labels - a logical flag. Should labels be printed?
	extend = 1.5
	u = x@parameter$u
	data = as.vector(x@data$exceedances)
	sorted = sort(data)
	shape = xi = x@fit$par.ests["xi"]
	scale = beta = x@fit$par.est["beta"]
	ypoints = ppoints(sorted)
	U = max(sorted)*extend
	z = qgpd(seq(0, 1, length = 1000), xi, u, beta)
	z = pmax(pmin(z, U), u)
	y = pgpd(z, xi, u, beta)
	# Labels:
	if (labels) {
		xlab = "Fu(x-u)"
		ylab = "x [log Scale]"
		main = "Excess Distribution"
	} else {
		xlab = ylab = main = ""
	}
	# Plot:
	plot(x = sorted, y = ypoints,
			xlim = range(u, U), ylim = range(ypoints, y, na.rm = TRUE),
			main = main, xlab = xlab, ylab = ylab,
			log = "x", axes = TRUE,
			col = "steelblue", pch = 19, ...)
	lines(z[y >= 0], y[y >= 0], col = "brown")
	# Addon:
	if (labels) {
		u = signif (x@parameter$u, 3)
		text = paste("u =", u)
		mtext(text, side = 3, adj = 0, cex = 0.7)
		grid()
	}
	mtext(paste("spd  : GPD Tail Fit"), side = 4, adj = 0, padj=0, col = "gray", cex = 0.7)
	# Return Value:
	invisible()
}

# ------------------------------------------------------------------------------
.plot.gpdfit.2<-function(x, labels = TRUE, ...)
{   # A function implemented by Diethelm Wuertz
	# and modified by Alexios Ghalanos
	# Description:
	#   Tail of Underlying Distribution
	# Arguments:
	#   x - an object of class fGPDFIT
	#   labels - a logical flag. Should labels be printed?
	extend = 1.5
	u = x@parameter$u
	data = as.vector(x@data$x)
	sorted = sort(data[data > u])
	prob = x@fit$prob
	shape = xi = x@fit$par.ests["xi"]
	beta = x@fit$par.ests["beta"]
	scale = beta * (1-prob)^xi
	location = u - (scale*((1 - prob)^(-xi)-1))/xi
	# Labels:
	if (labels) {
		xlab = "x [log scale]"
		ylab = "1-F(x) [log scale]"
		main = "Tail of Underlying Distribution"
	} else {
		xlab = ylab = main = ""
	}
	# Plot:
	U = max(data) * extend
	ypoints = ppoints(sorted)
	ypoints = (1 - prob) * (1 - ypoints)
	z = qgpd(seq(0, 1, length = 1000), xi, u, beta)
	z = pmax(pmin(z, U), u)
	y = pgpd(z, xi, u, beta)
	y = (1 - prob) * (1 - y)
	plot(x = sorted, y = ypoints,
			xlim = range(u, U), ylim = range(ypoints, y[y>0], na.rm = TRUE),
			main = main, xlab = xlab, ylab = ylab,
			log = "xy", axes = TRUE,
			col = "steelblue", pch = 19, ...)
	# Line:
	lines(z[y >= 0], y[y >= 0], col = "brown")
	if (labels) grid()
	mtext(paste("spd  : GPD Tail Fit"), side = 4, adj = 0, padj=0, col = "gray", cex = 0.7)

	# Return Value:
	invisible(list(x = sorted, y = ypoints))
}


# ------------------------------------------------------------------------------


.plot.gpdfit.3<-function(x, labels = TRUE, ...)
{   # A function implemented by Diethelm Wuertz
	# and modified by Alexios Ghalanos
	# Description:
	#   Scatterplot of GPD Residuals
	# Arguments:
	#   x - an object of class fGPDFIT
	#   labels - a logical flag. Should labels be printed?
	residuals = x@residuals

	# Labels:
	if (labels) {
		ylab = "Residuals"
		xlab = "Ordering"
		main = "Scatterplot of Residuals"
	} else {
		xlab = ylab = main = ""
	}

	# Plot:
	plot(residuals,
			main = main, ylab = ylab, xlab = xlab,
			col = "steelblue", pch = 19, ...)
	lines(lowess(1:length(residuals), residuals), col = "brown")
	if (labels) grid()
	mtext(paste("spd  : GPD Tail Fit"), side = 4, adj = 0, padj=0, col = "gray", cex = 0.7)
	# Return Value:
	invisible(list(x = 1:(length(residuals)), y = residuals))
}


# ------------------------------------------------------------------------------


.plot.gpdfit.4<-function(x, labels = TRUE, ...)
{   # A function implemented by Diethelm Wuertz
	# and modified by Alexios Ghalanos
	# Description:
	#   Quantile-Quantile Plot of GPD Residuals
	# Arguments:
	#   x - an object of class fGPDFIT
	#   labels - a logical flag. Should labels be printed?
	data = x@residuals
	sorted = sort(data)
	# Labels:
	if (labels) {
		xlab = "Ordered Data"
		ylab = "Exponential Quantiles"
		main = "QQ-Plot of Residuals"
	} else {
		xlab = ylab = main = ""
	}
	# Plot:
	y = qexp(ppoints(data))
	plot(x = sorted, y = y,
			main = main, xlab = xlab, ylab = ylab,
			col = "steelblue", pch = 19, ...)
	abline(lsfit(sorted, y), col = "brown")
	if (labels) grid()
	mtext(paste("spd  : GPD Tail Fit"), side = 4, adj = 0, padj=0, col = "gray", cex = 0.7)
	# Return Value:
	invisible(list(x = sorted, y = y))
}
