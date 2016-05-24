# LinearErrorsInVariablesClass.R
# Defines the class leiv for linear errors-in-variables objects
# Defines print and plot methods for the class
# Defines a partition integrate function
# Defines a median function
# Defines a probability interval function
# Author: David Leonard
# Date: 03 Jun 2012

# Date: 17 May 2010
# fix error when leiv is called with cor = 0

# Date: 21 May 2010
# fix error when leiv is called with zero covariance data
# fix special handling of singular cases in plot method

# Date: 26 May 2010
# replace beta cdf with equivalent F cdf

# Date: 28 Feb 2011
# stabilize probInt and p50 calculations

# Date: 20 Mar 2011
# remove priorType option

# Date: 26 Mar 2011
# improve efficiency of numerical integration by
# making direct calls to QUADPACK routines dqags and dqagi
# check normalization and add stop messages on failure

# Date: 30 Mar 2011
# improve checking for exceptional inputs

# Date: 09 Jun 2011
# add prior option to use another leiv object as the
# prior density

# Date: 03 Jun 2012
# replace direct calls to QUADPACK routines dqags and dqagi
# with calls to the wrapper function integrate

# median function
p50 <- function (p,interval,subdivisions=100,rel.tol=.Machine$double.eps^0.25,abs.tol=rel.tol) {
	# returns the median of the density p
	if (interval[1] > 0)
		prob <- function(x) integrate(p,x,interval[2],subdivisions=subdivisions,rel.tol=rel.tol,abs.tol=abs.tol)$value+integrate(p,interval[2],Inf,subdivisions=subdivisions,rel.tol=rel.tol,abs.tol=abs.tol)$value else
		prob <- function(x) integrate(p,-Inf,x,subdivisions=subdivisions,rel.tol=rel.tol,abs.tol=abs.tol)$value+integrate(p,interval[1],x,subdivisions=subdivisions,rel.tol=rel.tol,abs.tol=abs.tol)$value
	return(uniroot(function(x) prob(x)-0.5,interval,tol=abs.tol)$root)
}

# probability interval function
probInt <- function (p,interval,level,subdivisions=100,rel.tol=.Machine$double.eps^0.25,abs.tol=rel.tol) {
	# returns the shortest (100*level)% probability interval of the density p
	
	# find the mode
	opt <- optimize(p,interval,tol=abs.tol,maximum=TRUE)
	mode <- opt$maximum
	
	p12 <- function(interval) {
		if (interval[1] < 0 && interval[2] > 0) {
			partition <- c(interval[1],range(0,mode),interval[2])
			return(integrate(p,partition[1],partition[2],subdivisions=subdivisions,rel.tol=rel.tol,abs.tol=abs.tol)$value+integrate(p,partition[2],partition[3],subdivisions=subdivisions,rel.tol=rel.tol,abs.tol=abs.tol)$value+ integrate(p,partition[3],partition[4],subdivisions=subdivisions,rel.tol=rel.tol,abs.tol=abs.tol)$value)
		} else return(integrate(p,interval[1],mode,subdivisions=subdivisions,rel.tol=rel.tol,abs.tol=abs.tol)$value+integrate(p,mode,interval[2],subdivisions=subdivisions,rel.tol=rel.tol,abs.tol=abs.tol)$value)
	}
	
	# confirm valid interval
	while (p12(interval) < (1+level)/2) {
		interval <- mode+c(-1,1)*(interval[2]-interval[1])
		opt <- optimize(p,interval,tol=abs.tol,maximum=TRUE)
		mode <- opt$maximum
	}
	
	xlim <- function(y) {
		xlower <- optimize(function(x) (p(x)-y)^2,c(interval[1],mode),tol=abs.tol)$minimum
		xupper <- optimize(function(x) (p(x)-y)^2,c(mode,interval[2]),tol=abs.tol)$minimum
		return(c(xlower,xupper))
	}

	prob <- function(y) p12(xlim(y))	
	
	ymax <- opt$objective
	yOpt <- uniroot(function(y) prob(y)-level,c(0,ymax),f.lower=1-level,f.upper=-level,tol=abs.tol)$root
	
	return(xlim(yOpt))
}


# class definition

setClass("leiv",
	representation(
		slope="numeric",
		intercept="numeric",
		slopeInt="numeric",
		interceptInt="numeric",
		density="function",
		n="numeric",
		cor="numeric",
		sdRatio="numeric",
		xMean="numeric",
		yMean="numeric",
		call="call",
		probIntCalc="logical",
		level="numeric",
		x="numeric",
		y="numeric"
	)
)

# generating function

leiv <-
function(formula, data, subset, prior=NULL, n=NULL, cor=NULL, sdRatio=NULL, xMean=0, yMean=0, probIntCalc=FALSE, level=0.95, subdivisions=100, rel.tol=.Machine$double.eps^0.25, abs.tol=0.1*rel.tol, ... ) {
	cl <- match.call()
	if (is.null(n) || is.null(cor) || is.null(sdRatio)) {
		# interpret the call
		mf <- match.call(expand.dots = FALSE)
		m <- match(c("formula", "data", "subset"), names(mf), 0L)
		mf <- mf[c(1L, m)]
		mf$drop.unused.levels <- TRUE
		mf[[1L]] <- as.name("model.frame")
		mf <- eval(mf, parent.frame())
		y <- model.response(mf, "numeric")
		if (NCOL(y) != 1L)
			stop("only one y variable is supported")
		if (any(is.infinite(y)))
			stop("requires finite y data")
		mt <- attr(mf, "terms")
		attr(mt,"intercept") <- 0 # drop intercept from x
		x <- model.matrix(mt, mf)
		if (ncol(x) != 1L)
			stop("only one x variable is supported")
		x <- as.vector(x)
		if (any(is.infinite(x)))
			stop("requires finite x data")
		
		# sufficient statistics
		n <- length(x)
		if (n < 2L)
			stop("requires n >= 2 data points")
		Sxy <- cov(x,y)
		Sxx <- var(x)
		Syy <- var(y)
		if (Sxx > 0 && Syy > 0)	cor <- Sxy/sqrt(Sxx*Syy) else {
			if (Sxx == 0 && Syy == 0)
				stop("requires n > 2 distinct data points")
			if (Syy == 0)
				stop("singular data: linear with zero slope")
			if (Sxx == 0)
				stop("singular data: linear with infinite slope")
		}
		if (!(cor > -1 && cor < 1))
			warning(paste("singular data: cor =",cor))
		sdRatio <- sqrt(Syy/Sxx)
		xMean <- mean(x)
		yMean <- mean(y)
	} else {
		if (n < 2L)
			stop("requires n >= 2 data points")
		if (sdRatio < 0)
			stop("requires sdRatio > 0")
		if (sdRatio == 0)
			stop("singular input: linear with zero slope")
		if (is.infinite(sdRatio))
			stop("singular input: linear with infinite slope")
		if (n == 2 && cor != -1 && cor != 1)
			stop("requires |cor| = 1 if n = 2")
		if (!(cor > -1 && cor < 1))
			warning(paste("singular input: cor =",cor))
		x <- numeric()
		y <- numeric()
	}
	
	if (!is.null(prior) && class(prior)!="leiv")
		stop("prior must be NULL or a leiv object")

	# all of the following in terms of the dimensionless slope

	if (cor > -1 && cor < 1) {

		# intermediates
		v <- n-1
		vp <- v+1
		vm <- v-1
		s <- sqrt((1-cor^2)/v)

		I <- function(b,r) {
			# central integral of dimensionless, scalar arguments
			tLower <- -r/s
			tUpper <- (b-r)/s
			tIntegrand <- function(t) {
				F <- vm/vp*(v+t^2)/(tUpper+t-2*tLower)/(tUpper-t)
				return(dt(t,v)*pf(F,vp,vm))
			}
			return(integrate(tIntegrand,tLower,tUpper,subdivisions=subdivisions,rel.tol=rel.tol,abs.tol=abs.tol)$value)
		}

		# pseudo-likelihood of dimensionless b vector
		J <- function(b) sapply(b, function(b) {
			if (b == 0) return(0) else {
				bSign <- sign(b)
				bAbs <- bSign*b
				rbSign <- cor*bSign
				return(I(bAbs,rbSign)+I(1/bAbs,rbSign))
			}
		})
		
		# prior density
		if (is.null(prior)) {
			p0 <- function(b) dt(b,df=1) # cauchy
		} else {
			p0 <- function(b) prior@density(b*prior@sdRatio)*prior@sdRatio
		}

		# unnormalized posterior density
		p1 <- function(b) p0(b)*J(b)

		# normalized posterior density
		if (cor > 0) {
			bb <- c(cor,1/cor)
			bbb <- c(0,bb)
		} else {
			if (cor < 0) {
				bb <- c(1/cor,cor)
				bbb <- c(bb,0)
			} else {
				bb <- c(-1,1)
				bbb <- c(-1,0,1)
			}
		}
		error.action <- function(message) if (message != "OK") {
			if (message=="the input is invalid") msg <- "invalid tolerance parameters" else
			if (message=="the integral is probably divergent") msg <- "integral did not converge" else
			if (message=="roundoff error is detected in the extrapolation table") msg <- "algorithm did not converge" else
			if (message=="extremely bad integrand behaviour") msg <- "singular behavior detected" else
			if (message=="roundoff error was detected") msg <- "roundoff error detected" else
			msg <- "maximum number of subdivisions reached"
			stop(paste("posterior density normalization failed;",msg))
		}
		k1 <- integrate(p1,-Inf,bbb[1],subdivisions=subdivisions,rel.tol=rel.tol,abs.tol=abs.tol)
		error.action(k1$message)
		k2 <- integrate(p1,bbb[1],bbb[2],subdivisions=subdivisions,rel.tol=rel.tol,abs.tol=abs.tol)
		error.action(k2$message)
		k3 <- integrate(p1,bbb[2],bbb[3],subdivisions=subdivisions,rel.tol=rel.tol,abs.tol=abs.tol)
		error.action(k3$message)
		k4 <- integrate(p1,bbb[3],Inf,subdivisions=subdivisions,rel.tol=rel.tol,abs.tol=abs.tol)
		error.action(k4$message)
		k.value <- k1$value+k2$value+k3$value+k4$value
		if (k1$abs.err+k2$abs.err+k3$abs.err+k4$abs.err > k.value*rel.tol)
			stop("nonnormalized posterior density; try reducing abs.tol")
		p <- function(b) p1(b)/k.value

		# posterior median
		if (cor == 0) bMedian <- 0 else
			bMedian <- p50(p,bb,rel.tol=rel.tol,abs.tol=abs.tol)
	
		# probability interval
		if (probIntCalc) {
			blim <- bb-2/sqrt(n)*(bMedian-bb)
			bInt <- probInt(p,blim,level=level,rel.tol=rel.tol,abs.tol=abs.tol)
		} else {
			level <- numeric(0)
			bInt <- numeric(0)
		}

	} else {
		# exactly linear
		p <- function(b) ifelse(b == cor,1,0) # posterior density
		bMedian <- cor # posterior median
		if (probIntCalc) bInt <- c(cor,cor) else {
			level <- numeric(0)
			bInt <- numeric(0)
		}
	}
	
	# new leiv object with values in original units
	new("leiv",
		slope=bMedian*sdRatio,
		intercept=yMean-bMedian*sdRatio*xMean,
		slopeInt=bInt*sdRatio,
		interceptInt=if (xMean>0) (yMean-bInt*sdRatio*xMean)[c(2,1)] else yMean-bInt*sdRatio*xMean,
		density=function(b) p(b/sdRatio)/sdRatio,
		n=n,
		cor=cor,
		sdRatio=sdRatio,
		xMean=xMean,
		yMean=yMean,
		call=cl,
		probIntCalc=probIntCalc,
		level=level,
		x=x,
		y=y
	)
}

# print method

setMethod("print",
	signature(x="leiv"),
	function (x, digits = max(3, getOption("digits") - 3), ...) 
	{
    	cat("\nCall:\n", deparse(x@call), "\n", sep = "")
    	suffStats <- c(x@n,format(x@xMean,digits=digits),format(x@yMean,digits=digits),format(x@cor,digits=digits),format(x@sdRatio,digits=digits))
    	names(suffStats) <- c("sample size","mean x","mean y","sample cor", "sd ratio")
    	cat("\nSufficient statistics:\n")
    	print(suffStats,quote=FALSE)
    	cat("\n\nSlope estimate:",format(x@slope,digits=digits))
    	cat("\n\nIntercept estimate:",format(x@intercept,digits=digits))
		if (x@probIntCalc)
			cat("\n\nShortest",100*x@level,"% probability intervals:","\nslope (",toString(format(x@slopeInt,digits=digits)),")\nintercept (",toString(format(x@interceptInt,digits=digits)),")")
    	cat("\n\n")
    	invisible(x)
	}
)

# plot method

setMethod("plot",
	signature(x="leiv", y="missing"),
	function (x, plotType="density", xlim=NULL, ylim=NULL, xlab=NULL, ylab=NULL, col=NULL, lwd=NULL, ...) 
	{
		if (plotType == "scatter") {
			if (length(x@x) == 0L)
				stop("requires (x,y) data")
			if (is.null(xlab)) xlab <- "x"
			if (is.null(ylab)) ylab <- "y"
			if (is.null(col)) col <- "black"
			if (is.null(lwd)) lwd <- 1
			plot(x@x, x@y, xlim=xlim, ylim=ylim, xlab=xlab, ylab=ylab, ...)
			if (is.finite(x@slope)) abline(x@intercept, x@slope, col=col, lwd=lwd, ...) else
				abline(v=x@xMean, col=col, lwd=lwd, ...)
		} else {
			if (is.infinite(x@slope))
				stop("point mass at infinity")
			if (is.null(xlab)) xlab <- "slope"
			if (is.null(ylab)) ylab <- "density"
			if (is.null(col)) col <- "black"
			if (x@cor > -1 && x@cor < 1) {
				if (is.null(xlim) || is.null(ylim)) {
					if (x@cor == 0) {
						if (is.null(xlim)) xlim <- c(-1,1)-20/sqrt(x@n)*(x@slope-c(-1,1))
						if (is.null(ylim)) ylim <- c(0,1.2*optimize(x@density,xlim,maximum=TRUE)$objective)
					} else {
						if (x@cor > 0) bb <- c(x@cor,1/x@cor)*x@sdRatio else
							bb <- c(1/x@cor,x@cor)*x@sdRatio
						if (is.null(xlim)) xlim <- bb-20/sqrt(x@n)*(x@slope-bb)
						if (is.null(ylim)) ylim <- c(0,1.2*optimize(x@density,bb,maximum=TRUE)$objective)
					}
				}
				if (is.null(lwd)) lwd <- 2
				plot(x@density, xlim=xlim, ylim=ylim, xlab=xlab, ylab=ylab, col=col, lwd=lwd, ...)
			} else {
				if (is.null(xlim)) xlim <- x@slope+c(-1,1)
				if (is.null(ylim)) ylim <- c(0,1.2)
				if (is.null(lwd)) lwd <- 1
				plot(x=x@slope, y=1, xlim=xlim, ylim=ylim, xlab=xlab, ylab=ylab, col=col, lwd=lwd, ...)
			}
		}
   	invisible(x)
	}
)
