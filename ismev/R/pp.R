# This file contains the following functions:
# pp.fitrange  pp.fit  pp.diag  pp.pp  pp.qq
# ppf  ppq  ppp

"pp.fitrange"<-
function(data, umin, umax, npy = 365, nint = 10, show = FALSE, ...)
{
#
# produces estimates and 95% confidence intervals
# for point process model across range of thresholds
#
        m <- s <- up <- ul <- matrix(0, nrow = nint, ncol = 3)
	u <- seq(umin, umax, length = nint)
	for(i in 1:nint) {
		z <- pp.fit(data, u[i], npy, show = show, ...)
		m[i,  ] <- z$mle
		s[i,  ] <- z$se
		up[i,  ] <- z$mle + 1.96 * z$se
		ul[i,  ] <- z$mle - 1.96 * z$se
	}
	names <- c("Location", "Scale", "Shape")
	oldpar <- par(mfrow = c(1, 3))
	for(i in 1:3) {
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

"pp.fit"<-
function(xdat, threshold, npy = 365, ydat = NULL, mul = NULL, sigl = NULL, shl = NULL, mulink = identity, siglink = identity, shlink = identity, muinit = NULL, siginit = NULL, shinit = NULL, show = TRUE, method = "Nelder-Mead", maxit = 10000, ...)
{
	z <- list()
        npmu <- length(mul) + 1
	npsc <- length(sigl) + 1
	npsh <- length(shl) + 1
        n <- length(xdat)
	z$trans <- FALSE
	if(is.function(threshold)) 
            stop("`threshold' cannot be a function")
	u <- rep(threshold, length.out = n)
	if(length(unique(u)) > 1) z$trans <- TRUE
	uInd <- xdat > u
	lrate <- sum(uInd)/n
    	xdatu <- xdat[uInd]
	# xdatu <- xdat[xdat > u]
	# xind <- (1:n)[xdat > u]
	# u <- u[xind]
	in2 <- sqrt(6 * var(xdatu))/pi
	in1 <- mean(xdatu) - 0.57722 * in2
	if(is.null(shinit)) in3 <- 1e-8
	else in3 <- shinit
	in2 <- exp(log(in2) + in3*log(lrate))
	in1 <- threshold - (in2/in3)*(lrate^(-in3) - 1)
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
	init <- c(muinit, siginit, shinit)
	z$model <- list(mul, sigl, shl)
	z$link <- deparse(substitute(c(mulink, siglink, shlink)))
        z$threshold <- threshold
	z$npy <- npy
	z$nexc <- length(xdatu)
	z$data <- xdatu
        pp.lik <- function(a) {
	mu <- mulink(mumat %*% (a[1:npmu]))
	sc <- siglink(sigmat %*% (a[seq(npmu + 1, length = npsc)]))
	xi <- shlink(shmat %*% (a[seq(npmu + npsc + 1, length = npsh)]))
        if(any(sc^uInd <= 0)) return(10^6)
	if(min((1 + ((xi * (u - mu))/sc))^uInd) < 0) {
		l <- 10^6
	}
	else {
		y <- (xdat - mu)/sc
		y <- 1 + xi * y
		if(min(y) <= 0)
			l <- 10^6
		else l <- sum(uInd*log(sc)) + sum(uInd*log(y) * (1/xi + 1)) + n/npy * 
				mean((1 + (xi * (u - mu))/sc)^(-1/xi))
	}
	l
        }
	x <- optim(init, pp.lik, hessian = TRUE, method = method,
                   control = list(maxit = maxit, ...))
        mu <- mulink(mumat %*% (x$par[1:npmu]))
	sc <- siglink(sigmat %*% (x$par[seq(npmu + 1, length = npsc)]))
	xi <- shlink(shmat %*% (x$par[seq(npmu + npsc + 1, length = npsh)]))
	z$conv <- x$convergence
	z$nllh <- x$value
	z$vals <- cbind(mu, sc, xi, u)
	z$gpd <- apply(z$vals, 1, ppp, npy)
	if(z$trans) {
		z$data <- as.vector((1 + (xi[uInd] * (xdatu - u[uInd]))/z$gpd[2, uInd])^(-1/xi[uInd]))
	}
	z$mle <- x$par
        z$cov <- solve(x$hessian)
	z$se <- sqrt(diag(z$cov))
        if(show) {
	    if(z$trans)
		print(z[c(2, 3)])
	    if(length(z[[4]]) == 1)
		print(z[4])
	    print(z[c(5, 6, 8)])
	    if(!z$conv)
		print(z[c(9, 12, 14)])
        }
	class( z) <- "pp.fit"
        invisible(z)
}

"pp.diag"<-
function(z)
{
	n <- length(z$data)
	x <- (1:n)/(n + 1)
	if(z$trans) {
		oldpar <- par(mfrow = c(1, 2))
		plot(x, sort(z$data), xlab = "empirical", ylab = "model")
		abline(0, 1, col = 3)
		title("Residual Probability Plot")
		# plot( - log(1 - x),  - log(1 - sort(z$data)), ylab = "empirical", xlab = "model")
		plot( - log(1 - x),  sort(-log(z$data)), ylab = "empirical", xlab = "model")
		abline(0, 1, col = 3)
		title("Residual quantile Plot (Exptl. Scale)")
	}
	else {
		oldpar <- par(mfrow = c(1, 2), pty = "s")
		pp.pp(z$mle, z$threshold, z$npy, z$data)
		pp.qq(z$mle, z$threshold, z$npy, z$data)
	}
	par(oldpar)
	invisible()
}

"pp.pp"<-
function(a, u, npy, dat)
{
#
# function called by pp.diag
# produces probability plot
#
	y <- apply(as.matrix(sort(dat)), 1, ppf, a = a, u = u, npy = npy)
	plot((1:length(dat))/length(dat), y, xlab = "empirical", ylab = "model",
		main = "Probability plot")
	abline(0, 1, col = 4)
}

"pp.qq"<-
function(a, u, npy, dat)
{
#
# function called by pp.diag
# computes quantile plot
#
	y <- apply(as.matrix((length(dat):1/(length(dat) + 1))), 1, ppq, a = a, 
		u = u, npy = npy)
	plot(y, sort(dat), ylab = "empirical", xlab = "model", main = 
		"Quantile Plot")
	abline(0, 1, col = 4)
}

"ppf"<-
function(a, z, u, npy)
{
#
# ancillary function
# calculates distribution function in point process model
#
	b <- ppp(c(a, u), npy)
	1 - (1 + (b[3] * (z - u))/b[2])^(-1/b[3])
}

"ppq"<-
function(a, u, npy, p)
{
#
# ancillary routine
# finds quantiles in point process model
#
	b <- ppp(c(a, u), npy)
	u + (b[2] * (((p))^( - b[3]) - 1))/b[3]
}

"ppp"<-
function(a, npy)
{
	u <- a[4]
	la <- 1 - exp( - (1 + (a[3] * (u - a[1]))/a[2])^(-1/a[3])/npy)
	sc <- a[2] + a[3] * (u - a[1])
	xi <- a[3]
	c(la, sc, xi)
}

