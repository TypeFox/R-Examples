bwplot.numeric <- function(x, y, smax=13.49/IQR(x), normalize = FALSE, n.points = 1000, c.thresh = 2, ...)
{	
	n <- length(x)
	
	if( missing(y) ) {
		LINE = TRUE
		y <- rep(1, length(x))
		normalize = FALSE #Density estimates automatically normalize
	} else {
		LINE = FALSE
	}
	
	thresh <- c.thresh*sqrt(log(n, 10)/n)
	
	dft <- function(s) {
		dftval <- complex( real = sum(y * cos(s * x)), imaginary = sum(y * sin(s * x)) )
		if( normalize == TRUE ){
			return(Mod(dftval)/(n*abs(sum(y))))
		}
		else if (normalize == FALSE ){
			return(Mod(dftval)/n)
		}
	}

	svals <- seq(0, smax, length.out=n.points)
	ftvals <- sapply( svals, dft )
		
	plot( svals, ftvals, ylim=c(0, max(ftvals)), type="l", xlab ="s", ylab="f(s)")
	if( LINE & !is.null(c.thresh) ) abline(h=thresh, col="blue", lty=3)
}

bwplot.ts <- function(x, lag.max=NULL, c.thresh = 2, ...) {
	n <- length(x)
	thresh <- c.thresh*sqrt(log(n, 10)/n)
	
	ac <- acf(x, type="correlation", plot=FALSE, lag.max = lag.max)
	ac$acf <- abs(ac$acf)
	plot(ac, ci=0, ...)

	abline(h=thresh, col="blue", lty=3)
}

bwplot <- function(x, ...) UseMethod("bwplot")
