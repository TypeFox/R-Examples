"plot.abund" <-
function(x, n=x$n, lvert=TRUE, lvars=TRUE, lcol=2, llty=2, all=TRUE, dlab=c("cumsum", "% log(ind.)", "% non-zero"), dcol=c(1,2,4), dlty=c(par("lty"), par("lty"), par("lty")), dpos=c(1.5, 20), type="l", xlab="variables", ylab="abundance", main=paste("Abundance sorting for:",x$data, "with f =", round(x$f, 4)), ...) {
	# The following function actually draws the graph
	abund.graph <- function(X, N, Lvert, Lvars, Lcol, Llty, All, Dlab, Dcol, Dlty, Dpos, Type, Xlab, Ylab, Main, ...) {
		p <- length(X$p.log.ind)
		plot(X$p.log.ind, type="n", ylim=c(0, 100), xlab=Xlab, ylab=Ylab, main=Main, xaxs="i", xaxt="n", ...)
		axis(1, 1:p, labels=as.character(X$vr))
		# Do we plot all lines or not?
		if (All == FALSE) {
			lines(1:p, X$cumsum, col=Dcol[1], lty=Dlty[1], type=Type)
			# Since there is only one line, we don't need a legend!
		} else {
			lines(1:p, X$p.log.ind, col=Dcol[2], lty=Dlty[2], type=Type)
			lines(1:p, X$p.nonull, col=Dcol[3], lty=Dlty[3], type=Type)
			lines(1:p, X$cumsum, col=Dcol[1], lty=Dlty[1], type=Type)
			# Draw the legend
			legend(Dpos[1], Dpos[2], Dlab, col=Dcol, lty=Dlty)
		}
		if (is.null(N)==FALSE) { # We draw the lines
			# Verify N
			if (N < 0) N <- 0
			if (N > p) N <- p
			if (Lvert==TRUE)		# We draw a vertical line
				lines(c(N+0.5, N+0.5), c(-10,110), lty=Llty, col=Lcol)
			if (Lvars==TRUE) {		# We change colors of selected variables labels
				axis(1, 1:N, labels=as.character(X$vr[1:N]), col.axis=Lcol)
			}
		}
	}
	invisible(abund.graph(x, n, lvert, lvars, lcol, llty, all, dlab, dcol, dlty, dpos, type, xlab, ylab, main, ...))
}
