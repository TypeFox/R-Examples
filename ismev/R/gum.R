# This file contains the following functions:
# gum.fit  gum.diag  gum.rl  gum.df  gum.q  gum.dens

"gum.fit"<-
function(xdat, ydat = NULL, mul = NULL, sigl = NULL, mulink = identity, siglink = identity, muinit = NULL, siginit = NULL, show = TRUE, method = "Nelder-Mead", maxit = 10000, ...)
{
#
# finds mles etc for gumbel model
#
	z <- list()
        npmu <- length(mul) + 1
        npsc <- length(sigl) + 1
	z$trans <- FALSE
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
	z$model <- list(mul, sigl)
	z$link <- c(deparse(substitute(mulink)), deparse(substitute(siglink)))
	init <- c(muinit, siginit)
        gum.lik <- function(a) {
        # calculates neg log lik of gumbel model
	mu <- mulink(mumat %*% (a[1:npmu]))
	sc <- siglink(sigmat %*% (a[seq(npmu + 1, length = npsc)]))
        if(any(sc <= 0)) return(10^6)
	y <- (xdat - mu)/sc
        sum(log(sc)) + sum(y) + sum(exp( - y))
        }
	x <- optim(init, gum.lik, hessian = TRUE, method = method,
                   control = list(maxit = maxit, ...))
	z$conv <- x$convergence
        if(!z$conv) {
                mu <- mulink(mumat %*% (x$par[1:npmu]))
	        sc <- siglink(sigmat %*% (x$par[seq(npmu + 1, length = npsc)]))
	        z$nllh <- x$value
	        z$data <- xdat
	        if(z$trans) {
		        z$data <- as.vector((xdat - mu)/sc)
	        }
	        z$mle <- x$par
                z$cov <- solve(x$hessian)
	        z$se <- sqrt(diag(z$cov))
	        z$vals <- cbind(mu, sc)
        }
        if(show) {
	    if(z$trans)
		print(z[c(2, 3, 4)])
	    else print(z[4])
	    if(!z$conv)
                print(z[c(5, 7, 9)])
        }
	class( z) <- "gum.fit"
	invisible(z)
}

"gum.diag"<-
function(z)
{
#
# produces diagnostic plots for output of
# gum.fit stored in z
#
	z$mle <- c(z$mle, 0)
	n <- length(z$data)
	x <- (1:n)/(n + 1)
	if(z$trans) {
	        oldpar <- par(mfrow = c(1, 2))
	        plot(x, exp( - exp( - sort(z$data))), xlab = "empirical",
                     ylab = "model")
	       	abline(0, 1, col = 4)
	       	title("Residual Probability Plot")
	       	plot( - log( - log(x)), sort(z$data), xlab = 
	       		"empirical", ylab = "model")
	       	abline(0, 1, col = 4)
	       	title("Residual Quantile Plot (Gumbel Scale)")
	}
       	else {
       		oldpar <- par(mfrow = c(2, 2))
       		gev.pp(z$mle, z$data)
       		gev.qq(z$mle, z$data)
       		gum.rl(z$mle, z$cov, z$data)
       		gev.his(z$mle, z$data)
       	}
       	par(oldpar)
       	invisible()
}

"gum.rl"<-
function(a, mat, dat)
{
#
# function called by gum.diag
# produces return level curve and 95 % confidence intervals
# on usual scale for gumbel model
#
	eps <- 1e-006
	a1 <- a
	a2 <- a
	a1[1] <- a[1] + eps
	a2[2] <- a[2] + eps
	f <- c(seq(0.01, 0.09, by = 0.01), 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 
		0.8, 0.9, 0.95, 0.99, 0.995, 0.999)
	q <- gevq(a, 1 - f)
	d1 <- (gevq(a1, 1 - f) - q)/eps
	d2 <- (gevq(a2, 1 - f) - q)/eps
	d <- cbind(d1, d2)
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

"gum.df"<-
function(x, a, b)
{
#
# ancillary function calculates dist fnc of gumbel model
#
	exp( - exp( - (x - a)/b))
}

"gum.q"<-
function(x, a, b)
{
#
# ancillary routine
# calculates quantiles of gumbel distn
#
	a - b * log( - log(1 - x))
}

"gum.dens"<-
function(a, x)
{
#
# ancillary function calculates density for gumbel model
#
	y <- (x - a[1])/a[2]
	(exp( - y) * exp( - exp( - y)))/a[2]
}






