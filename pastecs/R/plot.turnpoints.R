"plot.turnpoints" <-
function(x, level=0.05, lhorz=TRUE, lcol=2, llty=2, type="l", xlab="data number", ylab=paste("I (bits), level = ", level*100, "%", sep=""), main=paste("Information (turning points) for:",x$data), ...) {
	# The next function actually draws the graph
	turnpoints.graph <- function(X, Level, Lhorz, Lcol, Llty, Type, Xlab, Ylab, Main, Sub, ...) {
		plot(X$tppos, X$info, type=Type, xlab=Xlab, ylab=Ylab, main=Main, ...)
		abline(h=-log(Level, base=2), lty=Llty, col=Lcol)
	}
	invisible(turnpoints.graph(x, level[1], lhorz, lcol, llty, type, xlab, ylab, main, ...))
}
