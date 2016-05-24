traceplot <-
function(x, 
	X=c("lambda","df"), 
#	labels=c("left", "right"),
	col = c("black", "red", "darkgreen", "blue","darkcyan","magenta", "brown","darkgray"), 
	pch = c(15:18, 7, 9, 12, 13),
	xlab, ylab="Coefficient", 
	xlim, ylim, ... ) {

	type <- X <- match.arg(X)
	if (type=="lambda") {
		X <- x$lambda
		if (missing(xlab)) xlab <- "Ridge constant"
		labels <- "left"
		}
	else {
		X <- x$df
		if (missing(xlab)) xlab <- "Degrees of freedom"
		labels <- "right"
	}
	coef <- coef(x)
	K <- nrow(coef)
	if (missing(xlim)) xlim <- range(X)
	if (missing(ylim)) ylim <- range(coef)

#	labels <- match.arg(labels)
	if (labels == "left") {
		xlim[1] <- xlim[1] - .1 * diff(xlim)
		labx <- X[1]
		laby <- coef[1,]
	}
	else {
		xlim[2] <- xlim[2] + .1 * diff(xlim)
		labx <- X[1]
		laby <- coef[1,]
	}

	matplot(X, coef, 	type="b", xlim=xlim, ylim=ylim, ylab=ylab, xlab=xlab, col=col, pch=pch, ...)
	abline(h=0, lty=3)
	if (type=="lambda") {
		abline(v=x$kHKB, col="gray", lty=2)
		text(x$kHKB, ylim[1], "HKB", pos=3)
		abline(v=x$kLW, col="gray", lty=2)
		text(x$kLW, ylim[1], "LW", pos=3)
	}
	vnames <- colnames(coef)
	text(labx, laby, colnames(coef), pos=c(2,4)[1+(labels=="right")])
}

