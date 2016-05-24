plot.TPCmsm <- function(x, plot.type="t", tr.choice, xlab, ylab, col, lty, xlim, ylim, conf.int=FALSE, ci.col, ci.lty, legend=TRUE, legend.pos, curvlab, legend.bty="n", ...) {
	if ( !inherits(x, "TPCmsm") ) stop("'x' must be of class 'TPCmsm'")
	if ( !( plot.type %in% c("t", "c") ) ) stop("Argument 'plot.type' must be one of 't' or 'c'")
	if ( missing(tr.choice) ) tr.choice <- dimnames(x$est)[[3]]
	lt <- length(tr.choice)
	if (sum( tr.choice %in% dimnames(x$est)[[3]] ) != lt) stop("Argument 'tr.choice' and possible transitions must match")
	if ( anyDuplicated(tr.choice) ) stop("Argument 'tr.choice' must be unique")
	lx <- length(x$x)
	itr <- match( tr.choice, dimnames(x$est)[[3]] )
	if (plot.type == "t") {
		if ( missing(xlim) ) xlim <- c(0, x$time[length(x$time)])
		if ( missing(ylim) ) ylim <- c(0, 1)
		if ( missing(xlab) ) xlab <- "Time"
		if ( missing(ylab) ) ylab <- "Transition probability"
		if ( missing(col) ) col <- rep(1, lt*lx)
		else if (length(col) < lt*lx) col <- col*rep(1, lt*lx)
		if ( missing(lty) ) lty <- seq_len(lt*lx)
		else if (length(lty) < lt*lx) lty <- lty*rep(1, lt*lx)
		plot(xlim, ylim, xlab=xlab, ylab=ylab, xlim=xlim, ylim=ylim, type="n", ...)
		for ( i in seq_len(lt) ) {
			for ( j in seq_len(lx) ) {
				lines(x=x$time, y=x$est[,which(x$covariate==x$x[j]),itr[i]], type="s", col=col[j+lx*(i-1)], lty=lty[j+lx*(i-1)], ...)
			}
		}
	} else if (plot.type == "c") {
		if ( missing(xlim) ) xlim <- c(x$covariate[1], x$covariate[length(x$covariate)])
		if ( missing(ylim) ) ylim <- c(0, 1)
		if ( missing(xlab) ) xlab <- "Covariate"
		if ( missing(ylab) ) ylab <- paste("P(", x$s, ", ", x$t, " | ", xlab, ")", sep="")
		if ( missing(col) ) col <- rep(1, lt)
		else if (length(col) < lt) col <- col*rep(1, lt)
		if ( missing(lty) ) lty <- seq_len(lt)
		else if (length(lty) < lt) lty <- lty*rep(1, lt)
		plot(xlim, ylim, xlab=xlab, ylab=ylab, xlim=xlim, ylim=ylim, type="n", ...)
		for ( i in seq_len(lt) ) {
			lines(x=x$covariate, y=x$est[dim(x$est)[1],,itr[i]], type="s", col=col[i], lty=lty[i], ...)
		}
	}
	if ( conf.int & !is.null(x$inf) & !is.null(x$sup) ) {
		if (plot.type == "t") {
			if ( missing(ci.col) ) ci.col <- col
			if (length(ci.col) < lt*lx) ci.col <- ci.col*rep(1, lt*lx)
			if ( missing(ci.lty) ) ci.lty <- 3
			if (length(ci.lty) < lt*lx) ci.lty <- ci.lty*rep(1, lt*lx)
			for ( i in seq_len(lt) ) {
				for ( j in seq_len(lx) ) {
					lines(x=x$time, y=x$inf[,which(x$covariate==x$x[j]),itr[i]], type="s", col=ci.col[j+lx*(i-1)], lty=ci.lty[j+lx*(i-1)], ...)
					lines(x=x$time, y=x$sup[,which(x$covariate==x$x[j]),itr[i]], type="s", col=ci.col[j+lx*(i-1)], lty=ci.lty[j+lx*(i-1)], ...)
				}
			}
		} else if (plot.type == "c") {
			if ( missing(ci.col) ) ci.col <- col
			if (length(ci.col) < lt) ci.col <- ci.col*rep(1, lt)
			if ( missing(ci.lty) ) ci.lty <- 3
			if (length(ci.lty) < lt) ci.lty <- ci.lty*rep(1, lt)
			for ( i in seq_len(lt) ) {
				lines(x=x$covariate, y=x$inf[dim(x$inf)[1],,itr[i]], type="s", col=ci.col[i], lty=ci.lty[i], ...)
				lines(x=x$covariate, y=x$sup[dim(x$sup)[1],,itr[i]], type="s", col=ci.col[i], lty=ci.lty[i], ...)
			}
		}
	}
	if (legend) {
		if ( missing(legend.pos) ) legend.pos <- "topleft"
		if ( missing(curvlab) ) {
			if (plot.type == "t" & lx > 1) {
				curvlab <- vector(mode="character", length=lt*lx)
				if (lt > 1) {
					for ( i in seq_len(lt) ) {
						for (j in seq_len(lx) ) {
							curvlab[j+lx*(i-1)] <- paste(tr.choice[i], "X =", x$x[j], sep=" ")
						}
					}
				} else {
					for ( j in seq_len(lx) ) {
						curvlab[j] <- paste("X =", x$x[j], sep=" ")
					}
				}
			} else {
				curvlab <- tr.choice
			}
		}
		addlegend(legend.pos, curvlab, col=col, lty=lty, legend.bty=legend.bty, ...)
	}
	invisible()
}
