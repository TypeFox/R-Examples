# Author: Robert J. Hijmans, r.hijmans@gmail.com
# Date: December 2009
# Version 0.1
# Licence GPL v3


if (!isGeneric("density")) {
	setGeneric("density", function(x, ...)
		standardGeneric("density"))
}	

setMethod('density', signature(x='DistModel'), 
	function(x, v=NULL, ...) {
		if (is.null(v)) {
			n <- ncol(x@presence)
			v <- 1:n
		} else {
			v <- round(v)
			v <- subset(v, v > 0)
			v <- subset(v, v <= ncol(x@presence))
			n <- length(v)
			if (n==0) {
				n <- ncol(x@presence)
				v <- 1:n				
			}
		}
		if (n > 1) {
			nc <- ceiling(sqrt(n))
			nr <- ceiling(n / nc)
			graphics::par(mfrow=c(nr, nc))
		}
		if (x@hasabsence) {
			for (i in 1:length(v)) {
				if (is.character(v[i])) {
					ab <- density(x@absence[,v[i]])
				} else {
					ab <- density(x@absence[,v[i]])
				}
				mt <- colnames(x@absence)[v[i]]
				pr <- density(x@presence[,i])
				yl = c(min(ab$y, pr$y), max(ab$y, pr$y))
				xl = c(min(ab$x, pr$x), max(ab$x, pr$x))
				plot(ab, main=mt, ylab='', xlab='', xlim=xl, ylim=yl, lwd=2, lty=2, col='blue', ...)
				lines(pr, col='red', lwd=2)
			} 
		} else {
			for (i in 1:length(v)) {
				if (is.character(v[i])) {
					plot(density(x@presence[,v[i]]), main=v[i], ylab='', xlab='', lwd=2, lty=2, col='red', ...)
				} else {
					plot(density(x@presence[,v[i]]), main=colnames(x@presence)[v[i]], ylab='', xlab='', lwd=2, lty=2, col='red', ...)
				}
			}
		}
	}
)
 
