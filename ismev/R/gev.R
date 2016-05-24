# This file contains the following functions:
# gev.fit  gev.diag  gev.pp  gev.qq  gev.rl  gev.his
# gevf  gevq  gev.dens  gev.profxi  gev.prof

"gev.fit"<-
function(xdat, ydat = NULL, mul = NULL, sigl = NULL, shl = NULL, mulink = identity, siglink = identity, shlink = identity, muinit = NULL, siginit = NULL, shinit = NULL, show = TRUE, method = "Nelder-Mead", maxit = 10000, ...)
{
#
# obtains mles etc for gev distn
#
	z <- list()
        npmu <- length(mul) + 1
        npsc <- length(sigl) + 1
        npsh <- length(shl) + 1
	z$trans <- FALSE	# if maximization fails, could try
# changing in1 and in2 which are 
# initial values for minimization routine
	in2 <- sqrt(6 * var(xdat))/pi
	in1 <- mean(xdat) - 0.57722 * in2
	if(is.null(mul)) {
		mumat <- as.matrix(rep(1, length(xdat)))
		if( is.null( muinit)) muinit <- in1
	}
	else {
		z$trans <- TRUE
		mumat <- cbind(rep(1, length(xdat)), ydat[, mul])
		if( is.null( muinit)) muinit <- c(in1, rep(0, length(mul)))
	}
	if(is.null(sigl)) {
		sigmat <- as.matrix(rep(1, length(xdat)))
		if( is.null( siginit)) siginit <- in2
	}
	else {
		z$trans <- TRUE
		sigmat <- cbind(rep(1, length(xdat)), ydat[, sigl])
		if( is.null( siginit)) siginit <- c(in2, rep(0, length(sigl)))
	}
	if(is.null(shl)) {
		shmat <- as.matrix(rep(1, length(xdat)))
		if( is.null( shinit)) shinit <- 0.1
	}
	else {
		z$trans <- TRUE
		shmat <- cbind(rep(1, length(xdat)), ydat[, shl])
		if( is.null( shinit)) shinit <- c(0.1, rep(0, length(shl)))
	}
	z$model <- list(mul, sigl, shl)
	z$link <- deparse(substitute(c(mulink, siglink, shlink)))
	init <- c(muinit, siginit, shinit)
        gev.lik <- function(a) {
        # computes neg log lik of gev model
        mu <- mulink(mumat %*% (a[1:npmu]))
        sc <- siglink(sigmat %*% (a[seq(npmu + 1, length = npsc)]))
	xi <- shlink(shmat %*% (a[seq(npmu + npsc + 1, length = npsh)]))
	y <- (xdat - mu)/sc
	y <- 1 + xi * y
	if(any(y <= 0) || any(sc <= 0)) return(10^6)
	sum(log(sc)) + sum(y^(-1/xi)) + sum(log(y) * (1/xi + 1))
        }
	x <- optim(init, gev.lik, hessian = TRUE, method = method,
                   control = list(maxit = maxit, ...))
	z$conv <- x$convergence
        mu <- mulink(mumat %*% (x$par[1:npmu]))
	sc <- siglink(sigmat %*% (x$par[seq(npmu + 1, length = npsc)]))
	xi <- shlink(shmat %*% (x$par[seq(npmu + npsc + 1, length = npsh)]))
	z$nllh <- x$value
	z$data <- xdat
	if(z$trans) {
		z$data <-  - log(as.vector((1 + (xi * (xdat - mu))/sc)^(
			-1/xi)))
	}
	z$mle <- x$par
        z$cov <- solve(x$hessian)
	z$se <- sqrt(diag(z$cov))
	z$vals <- cbind(mu, sc, xi)
        if(show) {
	    if(z$trans)
		print(z[c(2, 3, 4)])
	    else print(z[4])
	    if(!z$conv)
                print(z[c(5, 7, 9)])
	}
	class( z) <- "gev.fit"
	invisible(z)
}

"gev.diag"<-
function(z)
{
#
# produces diagnostic plots for output of
# gev.fit stored in z
#
	n <- length(z$data)
	x <- (1:n)/(n + 1)
	if(z$trans) {
      		oldpar <- par(mfrow = c(1, 2))
       		plot(x, exp( - exp( - sort(z$data))), xlab = 
       			"Empirical", ylab = "Model")
       		abline(0, 1, col = 4)
       		title("Residual Probability Plot")
       		plot( - log( - log(x)), sort(z$data), ylab = 
       			"Empirical", xlab = "Model")
       		abline(0, 1, col = 4)
       		title("Residual Quantile Plot (Gumbel Scale)")
       	}
       	else {
       		oldpar <- par(mfrow = c(2, 2))
       		gev.pp(z$mle, z$data)
       		gev.qq(z$mle, z$data)
       		gev.rl(z$mle, z$cov, z$data)
       		gev.his(z$mle, z$data)
       	}
       	par(oldpar)
       	invisible()
}

"gev.pp"<-
function(a, dat)
{
#
# sub-function for gev.diag
# produces probability plot
#
	plot((1:length(dat))/length(dat), gevf(a, sort(dat)), xlab = 
		"Empirical", ylab = "Model", main = "Probability Plot")
	abline(0, 1, col = 4)
}

"gev.qq"<-
function(a, dat)
{
#
# function called by gev.diag
# produces quantile plot
#
	plot(gevq(a, 1 - (1:length(dat)/(length(dat) + 1))), sort(dat), ylab = 
		"Empirical", xlab = "Model", main = "Quantile Plot")
	abline(0, 1, col = 4)
}

"gev.rl.gradient" <-
function (a, p) 
{
    # scale <- z$mle[2]
    # shape <- z$mle[3]
    scale <- a[2]
    shape <- a[3]
    if (shape < 0) 
        zero.p <- p == 0
    else zero.p <- logical(length(p))
    out <- matrix(NA, nrow = 3, ncol = length(p))
    out[1, ] <- 1
    if (any(zero.p)) {
        out[2, zero.p & !is.na(zero.p)] <- rep(-shape^(-1), sum(zero.p, 
            na.rm = TRUE))
        out[3, zero.p & !is.na(zero.p)] <- rep(scale * (shape^(-2)), 
            sum(zero.p, na.rm = TRUE))
    }
    if (any(!zero.p)) {
        yp <- -log(1 - p[!zero.p])
        out[2, !zero.p] <- -shape^(-1) * (1 - yp^(-shape))
        out[3, !zero.p] <- scale * (shape^(-2)) * (1 - yp^(-shape)) - 
            scale * shape^(-1) * yp^(-shape) * log(yp)
    }
    return(out)
}

"gev.rl"<-
function(a, mat, dat)
{
#
# function called by gev.diag
# produces return level curve and 95 % confidence intervals
# on usual scale
#
	eps <- 1e-006
	a1 <- a
	a2 <- a
	a3 <- a
	a1[1] <- a[1] + eps
	a2[2] <- a[2] + eps
	a3[3] <- a[3] + eps
	f <- c(seq(0.01, 0.09, by = 0.01), 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 
		0.8, 0.9, 0.95, 0.99, 0.995, 0.999)
	q <- gevq(a, 1 - f)
	# d1 <- (gevq(a1, 1 - f) - q)/eps
	# d2 <- (gevq(a2, 1 - f) - q)/eps
	# d3 <- (gevq(a3, 1 - f) - q)/eps
	# d <- cbind(d1, d2, d3)
	d <- t( gev.rl.gradient( a=a, p=1-f))
	v <- apply(d, 1, q.form, m = mat)
	plot(-1/log(f), q, log = "x", type = "n", xlim = c(0.1, 1000), ylim = c(
		min(dat, q), max(dat, q)), xlab = "Return Period", ylab = 
		"Return Level")
	title("Return Level Plot")
	lines(-1/log(f), q)
	lines(-1/log(f), q + 1.96 * sqrt(v), col = 4)
	lines(-1/log(f), q - 1.96 * sqrt(v), col = 4)
	points(-1/log((1:length(dat))/(length(dat) + 1)), sort(dat))
}

"gev.his"<-
function(a, dat)
{
#
# Plots histogram of data and fitted density
# for output of gev.fit stored in z
#
	h <- hist(dat, plot = FALSE)
	if(a[3] < 0) {
		x <- seq(min(h$breaks), min(max(h$breaks), (a[1] - a[2]/a[3] - 
			0.001)), length = 100)
	}
	else {
		x <- seq(max(min(h$breaks), (a[1] - a[2]/a[3] + 0.001)), max(h$
			breaks), length = 100)
	}
	y <- gev.dens(a, x)
	hist(dat, freq = FALSE, ylim = c(0, max(max(h$density),max(y))), xlab = "z", ylab = "f(z)", 
		main = "Density Plot")
	points(dat, rep(0, length(dat)))
	lines(x, y)
}

"gevf"<-
function(a, z)
{
#
# ancillary function calculates gev dist fnc
#
	if(a[3] != 0) exp( - (1 + (a[3] * (z - a[1]))/a[2])^(-1/a[3])) else 
			gum.df(z, a[1], a[2])
}

"gevq"<-
function(a, p)
{
	if(a[3] != 0)
		a[1] + (a[2] * (( - log(1 - p))^( - a[3]) - 1))/a[3]
	else gum.q(p, a[1], a[2])
}

"gev.dens"<-
function(a, z)
{
#
# evaluates gev density with parameters a at z
#
	if(a[3] != 0) (exp( - (1 + (a[3] * (z - a[1]))/a[2])^(-1/a[3])) * (1 + (
			a[3] * (z - a[1]))/a[2])^(-1/a[3] - 1))/a[2] else {
		gum.dens(c(a[1], a[2]), z)
	}
}

"gev.profxi"<-
function(z, xlow, xup, conf = 0.95, nint = 100)
{
#
# plots profile log-likelihood for shape parameter
# in gev model
#
	cat("If routine fails, try changing plotting interval", fill = TRUE)
	v <- numeric(nint)
	x <- seq(xup, xlow, length = nint)
	sol <- c(z$mle[1], z$mle[2])
        gev.plikxi <- function(a) {
        # computes profile neg log lik
        if (abs(xi) < 10^(-6)) {
                y <- (z$data - a[1])/a[2]
                if(a[2] <= 0) l <- 10^6
                else l <- length(y) * log(a[2]) + sum(exp(-y)) + sum(y)
        }
        else {
		y <- (z$data - a[1])/a[2]
		y <- 1 + xi * y
		if(a[2] <= 0 || any(y <= 0))
			l <- 10^6
		else l <- length(y) * log(a[2]) + sum(y^(-1/xi)) + sum(log(y
			)) * (1/xi + 1)
	}
	l
        }
	for(i in 1:nint) {
		xi <- x[i]
		opt <- optim(sol, gev.plikxi)
		sol <- opt$par ; v[i] <- opt$value
	}
	plot(x,  - v, type = "l", xlab = "Shape Parameter", ylab = 
		"Profile Log-likelihood")
	ma <-  - z$nllh
	abline(h = ma, col = 4)
	abline(h = ma - 0.5 * qchisq(conf, 1), col = 4)
	invisible()
}

"gev.prof"<-
function(z, m, xlow, xup, conf = 0.95, nint = 100)
{
#
# plots profile log likelihood for m 'year' return level
# in gev model
#
        if(m <= 1) stop("`m' must be greater than one")
	cat("If routine fails, try changing plotting interval", fill = TRUE)
	p <- 1/m
	v <- numeric(nint)
	x <- seq(xlow, xup, length = nint)
	sol <- c(z$mle[2], z$mle[3])
        gev.plik <- function(a) {
        # computes profile neg log lik
        if (abs(a[2]) < 10^(-6)) {
                mu <- xp + a[1] * log(-log(1 - p))
                y <- (z$data - mu)/a[1]
                if(is.infinite(mu) || a[1] <= 0) l <- 10^6
                else l <- length(y) * log(a[1]) + sum(exp(-y)) + sum(y)
        }
	else {
                mu <- xp - a[1]/a[2] * (( - log(1 - p))^( - a[2]) - 1)
		y <- (z$data - mu)/a[1]
		y <- 1 + a[2] * y
                if(is.infinite(mu) || a[1] <= 0 || any(y <= 0))
			l <- 10^6
		else l <- length(y) * log(a[1]) + sum(y^(-1/a[2])) + sum(log(
				y)) * (1/a[2] + 1)
	}
	l
        }
        for(i in 1:nint) {
                xp <- x[i]
		opt <- optim(sol, gev.plik)
		sol <- opt$par ; v[i] <- opt$value 
	}
	plot(x,  - v, type = "l", xlab = "Return Level", ylab = 
		" Profile Log-likelihood")
	ma <-  - z$nllh
	abline(h = ma, col = 4)
	abline(h = ma - 0.5 * qchisq(conf, 1), col = 4)
	invisible()
}







