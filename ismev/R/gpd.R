# This file contains the following functions:
# gpd.fitrange  gpd.fit  gpd.diag  gpd.pp  gpd.qq  gpd.rl
# gpd.his  gpdf  gpdq  gpdq2  gpd.dens  gpd.profxi  gpd.prof

"gpd.fitrange"<-
function(data, umin, umax, nint = 10, show = FALSE, ...)
{
#
# computes mle's in gpd model, adjusted for threshold, 
# over range of threshold choices.
#
	m <- s <- up <- ul <- matrix(0, nrow = nint, ncol = 2)
	u <- seq(umin, umax, length = nint)
	for(i in 1:nint) {
		z <- gpd.fit(data, u[i], show = show, ...)
		m[i,  ] <- z$mle
		m[i, 1] <- m[i, 1] - m[i, 2] * u[i]
		d <- matrix(c(1,  - u[i]), ncol = 1)
		v <- t(d) %*% z$cov %*% d
		s[i,  ] <- z$se
		s[i, 1] <- sqrt(v)
		up[i,  ] <- m[i,  ] + 1.96 * s[i,  ]
		ul[i,  ] <- m[i,  ] - 1.96 * s[i,  ]
	}
	names <- c("Modified Scale", "Shape")
	oldpar <- par(mfrow = c(2, 1))
	for(i in 1:2) {
		um <- max(up[, i])
		ud <- min(ul[, i])
		plot(u, m[, i], ylim = c(ud, um), xlab = "Threshold", ylab = 
			names[i], type = "b")
		for(j in 1:nint)
			lines(c(u[j], u[j]), c(ul[j, i], up[j, i]))
	}
        par(oldpar)
        invisible(list(thresholds=u, mle=m, se=s, ci.low=ul, ci.up=up))
}

"gpd.fit"<-
function(xdat, threshold, npy = 365, ydat = NULL, sigl = NULL, shl = NULL, siglink = identity, shlink = identity, siginit = NULL, shinit = NULL, show = TRUE, method = "Nelder-Mead", maxit = 10000, ...)
{
# 
# obtains mles etc for gpd model
#
	z <- list()
        npsc <- length(sigl) + 1
	npsh <- length(shl) + 1
        n <- length(xdat)
	z$trans <- FALSE
	if(is.function(threshold))
            stop("`threshold' cannot be a function")
	u <- rep(threshold, length.out = n)
        if(length(unique(u)) > 1) z$trans <- TRUE
	xdatu <- xdat[xdat > u]
	xind <- (1:n)[xdat > u]
	u <- u[xind]
	in2 <- sqrt(6 * var(xdatu))/pi
	in1 <- mean(xdatu, na.rm = TRUE) - 0.57722 * in2
	if(is.null(sigl)) {
		sigmat <- as.matrix(rep(1, length(xdatu)))
		if( is.null( siginit)) siginit <- in2
	}
	else {
		z$trans <- TRUE
		sigmat <- cbind(rep(1, length(xdatu)), ydat[xind, sigl])
		if( is.null( siginit)) siginit <- c(in2, rep(0, length(sigl)))
	}
	if(is.null(shl)) {
		shmat <- as.matrix(rep(1, length(xdatu)))
		if( is.null( shinit)) shinit <- 0.1
	}
	else {
		z$trans <- TRUE
		shmat <- cbind(rep(1, length(xdatu)), ydat[xind, shl])
		if( is.null( shinit)) shinit <- c(0.1, rep(0, length(shl)))
	}
	init <- c(siginit, shinit)
	z$model <- list(sigl, shl)
	z$link <- deparse(substitute(c(siglink, shlink)))
        z$threshold <- threshold
	z$nexc <- length(xdatu)
	z$data <- xdatu	
        gpd.lik <- function(a) {
        # calculates gpd neg log lik
	sc <- siglink(sigmat %*% (a[seq(1, length = npsc)]))
	xi <- shlink(shmat %*% (a[seq(npsc + 1, length = npsh)]))
	y <- (xdatu - u)/sc
	y <- 1 + xi * y
	if(min(sc) <= 0)
		l <- 10^6
	else {
		if(min(y) <= 0)
			l <- 10^6
		else {
			l <- sum(log(sc)) + sum(log(y) * (1/xi + 1))
		}
	}
	l
        }
        x <- optim(init, gpd.lik, hessian = TRUE, method = method,
                   control = list(maxit = maxit, ...))
	sc <- siglink(sigmat %*% (x$par[seq(1, length = npsc)]))
	xi <- shlink(shmat %*% (x$par[seq(npsc + 1, length = npsh)]))
	z$conv <- x$convergence
	z$nllh <- x$value
	z$vals <- cbind(sc, xi, u)
	if(z$trans) {
		z$data <-  - log(as.vector((1 + (xi * (xdatu - u))/sc)^(-1/xi))
			)
	}
	z$mle <- x$par
	z$rate <- length(xdatu)/n
        z$cov <- solve(x$hessian)
	z$se <- sqrt(diag(z$cov))
	z$n <- n
	z$npy <- npy
	z$xdata <- xdat
        if(show) {
	    if(z$trans)
		print(z[c(2, 3)])
	    if(length(z[[4]]) == 1)
		print(z[4])
	    print(z[c(5, 7)])
	    if(!z$conv)
		print(z[c(8, 10, 11, 13)])
        }
	class( z) <- "gpd.fit"
	invisible(z)
}

"gpd.diag"<-
function(z)
{
#
# produces diagnostic plots for gpd model
# estimated using gpd.fit with output stored in z
#
	n <- length(z$data)
	x <- (1:n)/(n + 1)
       	if(z$trans) {
       		oldpar <- par(mfrow = c(1, 2))
       		plot(x, 1 - exp( - sort(z$data)), xlab = "Empirical", 
       			ylab = "Model")
       		abline(0, 1, col = 4)
       		title("Residual Probability Plot")
       		plot( - log(1 - x), sort(z$data), ylab = "Empirical", 
       			xlab = "Model")
       		abline(0, 1, col = 4)
       		title("Residual Quantile Plot (Exptl. Scale)")
       	}
       	else {
       		oldpar <- par(mfrow = c(2, 2))
       		gpd.pp(z$mle, z$threshold, z$data)
       		gpd.qq(z$mle, z$threshold, z$data)
       		gpd.rl(z$mle, z$threshold, z$rate, z$n, z$npy, z$cov, z$
       			data, z$xdata)
       		gpd.his(z$mle, z$threshold, z$data)
       	}
        par(oldpar)
       	invisible()
}

"gpd.pp"<-
function(a, u, dat)
{
# 
# function called by gpd.diag
# produces probability plot for gpd model
#
	plot((1:length(dat))/length(dat), gpdf(a, u, sort(dat)), xlab = 
		"Empirical", ylab = "Model", main = "Probability Plot")
	abline(0, 1, col = 4)
}

"gpd.qq"<-
function(a, u, dat)
{
#
# function called by gpd.diag
# produces quantile plot for gpd model
#
	plot(gpdq(a, u, 1 - (1:length(dat)/(length(dat) + 1))), sort(dat), ylab
		 = "Empirical", xlab = "Model", main = "Quantile Plot")
	abline(0, 1, col = 4)
}

"gpd.rl.gradient" <-
function (a, m) 
{
    # 'a' is a numeric vector of length 4 containing Pr(X>u) and the 
    # estimated scale and shape parameters.
    # a <- c(z$rate, z$mle)
    if (length(a) > 3) 
        stop("gpd.rl.gradient: Covariates not allowed.")
    out <- matrix(0, nrow = 3, ncol = length(m))
    if (a[3] == 0) {
        out[1, ] <- a[2]/a[1]
        out[2, ] <- log(m * a[1])
    }
    else {
        out[1, ] <- a[2] * m^(a[3]) * a[1]^(a[3] - 1)
        out[2, ] <- a[3]^(-1) * ((m * a[1])^(a[3]) - 1)
        out[3, ] <- -a[2] * a[3]^(-2) * ((m * a[1])^(a[3]) - 
            1) + a[2] * a[3]^(-1) * (m * a[1])^(a[3]) * log(m * 
            a[1])
    }
    return(out)
}

"gpd.rl"<-
function(a, u, la, n, npy, mat, dat, xdat)
{
#
# function called by gpd.diag
# produces return level curve and 95% confidence intervals
# for fitted gpd model
	a <- c(la, a)
	eps <- 1e-006
	a1 <- a
	a2 <- a
	a3 <- a
	a1[1] <- a[1] + eps
	a2[2] <- a[2] + eps
	a3[3] <- a[3] + eps
	jj <- seq(-1, 3.75 + log10(npy), by = 0.1)
	m <- c(1/la, 10^jj)
	q <- gpdq2(a[2:3], u, la, m)
	# d1 <- (gpdq2(a1[2:3], u, la, m) - q)/eps
	# d2 <- (gpdq2(a2[2:3], u, la, m) - q)/eps
	# d3 <- (gpdq2(a3[2:3], u, la, m) - q)/eps
	# d <- cbind(d1, d2, d3)
	d <- t( gpd.rl.gradient( a=a, m=m))
	mat <- matrix(c((la * (1 - la))/n, 0, 0, 0, mat[1, 1], mat[1, 2], 0, 
		mat[2, 1], mat[2, 2]), ncol = 3)
	v <- apply(d, 1, q.form, m = mat)
	plot(m/npy, q, log = "x", type = "n", xlim = c(0.1, max(m)/npy), ylim
		 = c(u, max(xdat, q[q > u - 1] + 1.96 * sqrt(v)[q > u - 1])), 
		xlab = "Return period (years)", ylab = "Return level", main = 
		"Return Level Plot")
	lines(m[q > u - 1]/npy, q[q > u - 1])
	lines(m[q > u - 1]/npy, q[q > u - 1] + 1.96 * sqrt(v)[q > u - 1], col
		 = 4)
	lines(m[q > u - 1]/npy, q[q > u - 1] - 1.96 * sqrt(v)[q > u - 1], col
		 = 4)
	nl <- n - length(dat) + 1
	sdat <- sort(xdat)
	points((1/(1 - (1:n)/(n + 1))/npy)[sdat > u], sdat[sdat > u])	
	#	points(1/(1 - (1:n)/(n + 1))/npy, 
#		sort(xdat))
#	abline(h = u, col = 3)
}

"gpd.his"<-
function(a, u, dat)
{
#
# function called by gpd.diag
# produces histogram and density plot
#
	h <- hist(dat, plot = FALSE)
	x <- seq(u, max(h$breaks), length = 100)
	y <- gpd.dens(a, u, x)
	hist(dat, freq = FALSE, ylim = c(0, max(max(h$density), max(y))), xlab = "x", ylab = "f(x)", 
		main = "Density Plot")
	lines(x, y, col = 4)
}

"gpdf"<-
function(a, u, z)
{
#
# ancillary function
# calculates gpd distribution function
#
	1 - (1 + (a[2] * (z - u))/a[1])^(-1/a[2])
}

"gpdq"<-
function(a, u, p)
u + (a[1] * (p^( - a[2])	#
# ancillary function
# computes gpd quantiles
#
 - 1))/a[2]

"gpdq2"<-
function(a, u, la, m)
{
#
# ancillary function
# calculates quantiles of gpd model
#
	u + (a[1] * ((m * la)^(a[2]) - 1))/a[2]
}

"gpd.dens"<-
function(a, u, z)
{
#
# ancillary function computes gpd density
#
	(1 + (a[2] * (z - u))/a[1])^(-1/a[2] - 1)/a[1]
}

"gpd.profxi"<-
function(z, xlow, xup, conf = 0.95, nint = 100)
{
#
# plots profile log likelihood for shape parameter
# in gpd model
#
	cat("If routine fails, try changing plotting interval", fill = TRUE)
	xdat <- z$data ; u <- z$threshold
	v <- numeric(nint)
	x <- seq(xup, xlow, length = nint)
	sol <- z$mle[1]
        gpd.plikxi <- function(a) {
        # calculates profile log lik
	if(abs(xi) < 10^(-4)) l <- length(xdat) * log(a) + sum(xdat - u)/a
		 else {
		y <- (xdat - u)/a
		y <- 1 + xi * y
                if(any(y <= 0) || a <= 0)
			l <- 10^6
		else l <- length(xdat) * log(a) + sum(log(y)) * (1/xi + 1)
	}
	l
        }
	for(i in 1:nint) {
		xi <- x[i]
		opt <- optim(sol, gpd.plikxi, method = "BFGS")
		sol <- opt$par ; v[i] <- opt$value 
	}
	plot(x,  - v, type = "l", xlab = "Shape Parameter", ylab = 
		"Profile Log-likelihood")
	ma <-  - z$nllh
	abline(h = ma, lty = 1)
	abline(h = ma - 0.5 * qchisq(conf, 1), lty = 1)
	invisible()
}



"gpd.prof"<-
function(z, m, xlow, xup, npy = 365, conf = 0.95, nint = 100)
{
#
# plots profile log-likelihood for m-year return level
# in gpd model
#
	cat("If routine fails, try changing plotting interval", fill = TRUE)
        xdat <- z$data ; u <- z$threshold ; la <- z$rate
	v <- numeric(nint)
	x <- seq(xlow, xup, length = nint)
        m <- m * npy
	sol <- z$mle[2]
        gpd.plik <- function(a) {
        # calculates profile neg log lik
        if(m != Inf) sc <- (a * (xp - u))/((m * la)^a - 1) else sc <- (u - xp)/
			a
	if(abs(a) < 10^(-4))
		l <- length(xdat) * log(sc) + sum(xdat - u)/sc
	else {
		y <- (xdat - u)/sc
		y <- 1 + a * y
                if(any(y <= 0) || sc <= 0)
			l <- 10^6
		else l <- length(xdat) * log(sc) + sum(log(y)) * (1/a + 1)
	}
	l
        }
	for(i in 1:nint) {
		xp <- x[i]
		opt <- optim(sol, gpd.plik, method = "BFGS")
		sol <- opt$par ; v[i] <- opt$value
	}
	plot(x,  - v, type = "l", xlab = "Return Level", ylab = 
		"Profile Log-likelihood")
	ma <-  - z$nllh
	abline(h = ma)
	abline(h = ma - 0.5 * qchisq(conf, 1))
	invisible()
}





