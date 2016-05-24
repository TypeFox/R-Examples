features <- function(x, y, smoother=c("glkerns", "smooth.spline"), fits.return=TRUE, control=list( ), ...) 
{

ctrl <- list(npts=max(100, length(y)), c.outlier=3, decim.out=2)

namc <- names(control)
    if (!all(namc %in% names(ctrl))) 
        stop("unknown names in control: ", namc[!(namc %in% names(ctrl))])
    ctrl[namc] <- control
    npts <- ctrl$npts
    plot.it <- ctrl$plot.it
    c.outlier <- ctrl$c.outlier
    decim.out <- ctrl$decim.out

trapezoid <- function(x,y) sum(diff(x)*(y[-1] + y[-length(y)]))/2

smoother <- match.arg(smoother, c("glkerns", "smooth.spline"))

if (smoother == "glkerns") {
	deriv1 <- function(z, ...) glkerns(x, y, deriv=1, x.out=z, hetero=TRUE, ...)$est
	deriv2 <- function(z, ...) glkerns(x, y, deriv=2, x.out=z, hetero=TRUE, ...)$est
}

if (smoother == "smooth.spline") {
	deriv1 <- function(z, ...) predict(fit, deriv=1, x=z)$y
	deriv2 <- function(z, ...) predict(fit, deriv=2, x=z)$y
}

n <- length(x)
cp <- cv <- ol <- NA


x.out <- seq(min(x), max(x), length=max(npts, length(x)))

if (smoother == "glkerns") {
fit <- glkerns(x, y, deriv=0, x.out=x, hetero=FALSE, ...)
fit0 <- glkerns(x, y, deriv=0, x.out=x.out, hetero=TRUE, ...)$est
fit1 <- glkerns(x, y, deriv=1, x.out=x.out, hetero=TRUE, ...)$est
fit2 <- glkerns(x, y, deriv=2, x.out=x.out, hetero=TRUE, ...)$est
resid <- y - fit$est
resid.scaled <- abs(scale(resid))
if (fits.return) {
	fits <- fit$est
	fits1 <- glkerns(x, y, deriv=1, x.out=x, hetero=TRUE, ...)$est
	fits2 <- glkerns(x, y, deriv=2, x.out=x, hetero=TRUE, ...)$est
	}
}

if (smoother == "smooth.spline") {
fit <- smooth.spline(x, y, cv=FALSE, ...)
fit0 <- predict(fit, deriv=0, x=x.out)$y
fit1 <- predict(fit, deriv=1, x=x.out)$y
fit2 <- predict(fit, deriv=2, x=x.out)$y
resid <- y - predict(fit)$y
resid.scaled <- abs(scale(resid))
if (fits.return) {
	fits <- predict(fit)$y
	fits1 <- predict(fit, deriv=1)$y
	fits2 <- predict(fit, deriv=2)$y
	}
}

fmean <- trapezoid(x.out, fit0) / diff(range(x.out))
fstar <- sqrt(trapezoid(x.out, fit0^2) / diff(range(x.out)))
fsd <- sqrt(fstar^2 - fmean^2)
noise <- attr(resid.scaled, "scaled:scale")
snr <- fsd/noise
fwiggle <- sqrt(trapezoid(x.out, fit2^2) / diff(range(x.out)))

d0 <- (fit1[-1] * fit1[-npts]) < 0
ncrt <- sum(d0)
if (ncrt == 0) crtpts <- curv <- NULL  else {
	crtpts <- curv <- rep(NA, ncrt)
	ind0 <- (1:npts)[d0]
	for (i in 1:ncrt) {
	temp <- try(uniroot(interval=c(x.out[ind0[i]], x.out[1 + ind0[i]]), f=deriv1), silent=TRUE)
	if (class(temp) != "try-error") {
		crtpts[i] <- temp$root
		curv[i] <- deriv2(temp$root)
		}
	}
    }

rm(fit)
outl <- resid.scaled > c.outlier
if (sum(outl) == 0) outliers <- NULL else outliers <- x[outl] 

if (!is.null(crtpts) ) {
	cp <- crtpts
	cv <- curv
	}

if (!is.null(outliers) ) ol <- outliers

f <- as.numeric(c(fmean, as.numeric(range(fit0)), fsd, noise, snr, as.numeric(range(fit1)), fwiggle, length(crtpts)))

names(f) <- c("fmean", "fmin", "fmax", "fsd", "noise", "snr", "d1min", "d1max", "fwiggle", "ncpts")

cp <- round(cp, decim.out)
cv <- round(cv, decim.out)
ol <- round(ol, decim.out)

ret.obj <- list(f=round(f, decim.out), cpts=cp, curvature=cv, outliers=ol)
if (fits.return) attr(ret.obj, "fits") <- list(x=x, y=y, fn=fits, d1=fits1, d2=fits2)
class(ret.obj) <- "features"
return(ret.obj)
}

##############
plot.features <- function(x,  ...) {
if (class(x) != "features") stop("Input must be an object of class `features' ")
fits <- attr(x, "fits")
old.par <- par(no.readonly = TRUE)
on.exit(par(old.par)) 
par(mfrow=c(2,2))
plot(fits$x, fits$y, xlab="x", ylab="smoothed function", ...)
lines(fits$x, fits$fn, ...)
plot(fits$x, fits$y, xlab="x", ylab="smoothed function", ...)
lines(fits$x, fits$fn, ...)
plot(fits$x, fits$d1, type="l", xlab="x", ylab="First Derivative", ...)
abline(h=0, lty=3)
plot(fits$x, fits$d2, type="l", xlab="x", ylab="Second Derivative", ...)
abline(h=0, lty=3)
}


fget <- function(x) { 
if (class(x) != "features") stop("Input must be an object of class `features' ") else UseMethod("fget") 
}

fget.features <- function(x) {
if (class(x) != "features") stop("Input must be an object of class `features' ")
list(f=x$f, crit.pts=x$cpts, curvature=x$curvature, outlier=x$outlier)
}




