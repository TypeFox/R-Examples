plot.TPmsm <- function(x, tr.choice, xlab="Time", ylab="Transition probability", col, lty, xlim, ylim, conf.int=FALSE, ci.col, ci.lty, legend=TRUE, legend.pos, curvlab, legend.bty="n", ...) {
	if ( !inherits(x, "TPmsm") ) stop("'x' must be of class 'TPmsm'")
	if ( missing(tr.choice) ) tr.choice <- colnames(x$est)
	lt <- length(tr.choice)
	if (sum( tr.choice %in% colnames(x$est) ) != lt) stop("Argument 'tr.choice' and possible transitions must match")
	if ( anyDuplicated(tr.choice) ) stop("Argument 'tr.choice' must be unique")
	if ( missing (xlim) ) xlim <- c( 0, max(x$time) )
	if ( missing(ylim) ) ylim <- c(0, 1)
	if ( missing(col) ) col <- rep(1, lt)
	else if (length(col) < lt) col <- col*rep(1, lt)
	if ( missing(lty) ) lty <- seq_len(lt)
	else if (length(lty) < lt) lty <- lty*rep(1, lt)
	itr <- match( tr.choice, colnames(x$est) )
	plot(xlim, ylim, xlab=xlab, ylab=ylab, xlim=xlim, ylim=ylim, type="n", ...)
	for ( i in seq_len(lt) ) {
		lines(x=x$time, y=x$est[,itr[i]], type="s", col=col[i], lty=lty[i], ...)
	}
	if ( conf.int & !is.null(x$inf) & !is.null(x$sup) ) {
		if ( missing(ci.col) ) ci.col <- col
		if (length(ci.col) < lt) ci.col <- ci.col*rep(1, lt)
		if ( missing(ci.lty) ) ci.lty <- 3
		if (length(ci.lty) < lt) ci.lty <- ci.lty*rep(1, lt)
		for ( i in seq_len(lt) ) {
			lines(x=x$time, y=x$inf[,itr[i]], type="s", col=ci.col[i], lty=ci.lty[i], ...)
			lines(x=x$time, y=x$sup[,itr[i]], type="s", col=ci.col[i], lty=ci.lty[i], ...)
		}
	}
	if (legend) {
		if ( missing(legend.pos) ) legend.pos <- "topleft"
		if ( missing(curvlab) ) curvlab <- tr.choice
		addlegend(legend.pos, curvlab, col=col, lty=lty, legend.bty=legend.bty, ...)
	}
	invisible()
}
