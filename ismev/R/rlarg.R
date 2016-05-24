# This file contains the following functions:
# rlarg.fit  rlarg.diag  rlarg.pp  rlarg.qq
# rlargf  rlargq  rlargq2

"rlarg.fit"<-
function(xdat, r = dim(xdat)[2], ydat = NULL, mul = NULL, sigl = NULL, shl = NULL, mulink = identity, siglink = identity, shlink = identity, muinit = NULL, siginit = NULL, shinit = NULL, show = TRUE, method = "Nelder-Mead", maxit = 10000, ...)
{
#
# calculates mles etc for rlargest order statistic model
#
	z <- list()
        npmu <- length(mul) + 1
        npsc <- length(sigl) + 1
        npsh <- length(shl) + 1
        z$trans <- FALSE
	in2 <- sqrt(6 * var(xdat[, 1]))/pi
	in1 <- mean(xdat[, 1]) - 0.57722 * in2
	if(is.null(mul)) {
		mumat <- as.matrix(rep(1, dim(xdat)[1]))
		if( is.null( muinit)) muinit <- in1
	}
	else {
		z$trans <- TRUE
		mumat <- cbind(rep(1, dim(xdat)[1]), ydat[, mul])
		if( is.null( muinit)) muinit <- c(in1, rep(0, length(mul)))
	}
	if(is.null(sigl)) {
		sigmat <- as.matrix(rep(1, dim(xdat)[1]))
		if( is.null( siginit)) siginit <- in2
	}
	else {
		z$trans <- TRUE
		sigmat <- cbind(rep(1, dim(xdat)[1]), ydat[, sigl])
		if( is.null( siginit)) siginit <- c(in2, rep(0, length(sigl)))
	}
	if(is.null(shl)) {
		shmat <- as.matrix(rep(1, dim(xdat)[1]))
		if( is.null( shinit)) shinit <- 0.1
	}
	else {
		z$trans <- TRUE
		shmat <- cbind(rep(1, dim(xdat)[1]), ydat[, shl])
		if( is.null( shinit)) shinit <- c(0.1, rep(0, length(shl)))
	}
        xdatu <- xdat[, 1:r, drop = FALSE]
        init <- c(muinit, siginit, shinit)
	z$model <- list(mul, sigl, shl)
	z$link <- deparse(substitute(c(mulink, siglink, shlink)))
        u <- apply(xdatu, 1, min, na.rm = TRUE)
        rlarg.lik <- function(a) {
        # calculates neg log lik
	mu <- mulink(drop(mumat %*% (a[1:npmu])))
	sc <- siglink(drop(sigmat %*% (a[seq(npmu + 1, length = npsc)])))
	xi <- shlink(drop(shmat %*% (a[seq(npmu + npsc + 1, length = npsh)])))
        if(any(sc <= 0)) return(10^6)
        y <- 1 + xi * (xdatu - mu)/sc
	if(min(y, na.rm = TRUE) <= 0)
		l <- 10^6
	else {
                y <- (1/xi+1) * log(y) + log(sc)
                y <- rowSums(y, na.rm = TRUE)
                l <- sum((1 + xi * (u - mu)/sc)^(-1/xi) + y)
        }
	l
        }
	x <- optim(init, rlarg.lik, hessian = TRUE, method = method,
                   control = list(maxit = maxit, ...))
        mu <- mulink(drop(mumat %*% (x$par[1:npmu])))
	sc <- siglink(drop(sigmat %*% (x$par[seq(npmu + 1, length = npsc)])))
	xi <- shlink(drop(shmat %*% (x$par[seq(npmu + npsc + 1, length = npsh)])))
	z$conv <- x$convergence
	z$nllh <- x$value
	z$data <- xdat
	if(z$trans) {
		for(i in 1:r)
			z$data[, i] <-  - log((1 + (as.vector(xi) * (xdat[, i] - 
				as.vector(mu)))/as.vector(sc))^(-1/as.vector(xi
				)))
	}
	z$mle <- x$par
        z$cov <- solve(x$hessian)
	z$se <- sqrt(diag(z$cov))
	z$vals <- cbind(mu, sc, xi)
	z$r <- r
        if(show) {
	    if(z$trans)
		print(z[c(2, 3)])
	    print(z[4])
	    if(!z$conv)
		print(z[c(5, 7, 9)])
        }
	class( z) <- "rlarg.fit"
	invisible(z)
}

"rlarg.diag"<-
function(z, n = z$r)
{
#
# takes output from rlarg.fit
# produces probability and quantile plots for
# each order statistic
#
	z2 <- z
	z2$data <- z$data[, 1]
        oldpar <- par(ask = TRUE, mfcol = c(2, 2))
	if(z$trans) {
		for(i in 1:n) {
			rlarg.pp(c(0, 1, 0), z$data[, 1:z$r], i)
			rlarg.qq(c(0, 1, 0), z$data[, 1:z$r], i)
		}
	}
	else {
		gev.diag(z2)
		for(i in 1:n) {
			rlarg.pp(z$mle, z$data, i)
			rlarg.qq(z$mle, z$data, i)
		}
	}
	par(oldpar)
	invisible()
}

"rlarg.pp"<-
function(a, dat, k)
{
#
# ancillary function
# calculates probability plot in r largest model
#
	da <- dat[!is.na(dat[, k]), k]
	plot((1:length(da))/length(da), rlargf(a, sort(da), k), xlab = "", ylab
		 = "")
	title(paste("k=", k, sep = ""), cex = 0.7)
	abline(0, 1, col = 4)
}

"rlarg.qq"<-
function(a, dat, k)
{
#
# ancillary function
# calculates quantile plot in r largest model
#
	da <- dat[!is.na(dat[, k]), k]
	plot(rlargq(a, 1 - (1:length(da)/(length(da) + 1)), k, da), sort(da), 
		xlab = "", ylab = "")
	title(paste("k=", k, sep = ""), cex = 0.7)
	abline(0, 1, col = 4)
}

"rlargf"<-
function(a, z, k)
{
#
# ancillary function
# calculates dist fnc in r largest model
#
	eps <- 10^(-6)
	res <- NULL
	if(abs(a[3]) < eps)
		tau <- exp( - (z - a[1])/a[2])
	else tau <- (1 + (a[3] * (z - a[1]))/a[2])^(-1/a[3])
	for(i in 1:length(tau)) {
		if(is.na(tau[i]))
			res[i] <- 1
		else res[i] <- exp( - tau[i]) * sum(tau[i]^(0:(k - 1))/gamma(1:(
				k)))
	}
	res
}

"rlargq"<-
function(a, p, k, dat)
{
#
# ancillary routine 
# for finding quantiles in r largest model
	res <- NULL
	for(i in 1:length(p)) {
		inter <- c(min(dat) - 1, max(dat) + 1)
		res[i] <- uniroot(rlargq2, inter, a = a, kk = k, p = p[i])$root
	}
	res
}

"rlargq2"<-
function(x, a, kk, p)
{
#
# ancillary routine
# for finding quantiles in r largest model
#
	res <- rlargf(a, x, kk) - (1 - p)
	res
}



