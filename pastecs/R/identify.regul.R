"identify.regul" <-
function(x, series=1, col=3, label="#", ...) {
	labels <- rep(label, length.out=length(x$xini))
	i <- series
	if (i > ncol(x$yini))
		stop("This series does not exist")
	n <- identify(x$xini, x$yini[,i], labels=labels, col=col, ...)
	n.vec <- rep(0, length.out=length(x$xini))
	n.vec[n] <- 1
	n.vec
}
