plot.cdfDist <-
function(x, xlab='x', ylab='P(X < x)',
	rlab='Cumulative measure', col=c('#225588', '#88AA55', '#DD88AA'),
	lty=1:3, lwd=1+1:3/3, type='l', axes=TRUE, plotMu=TRUE, sigdig=2,
	leg=c("Observed", "Expected", "Cumulative measure"), ylim=0:1, ...){
	plot(x$x, x$F1, type='l', col=col[1], axes=FALSE, xlab=xlab,
			ylab=ylab, ylim=ylim, ...)
	lines(x$x, x$F2, col=col[2], lwd=lwd[2], lty=lty[2])
	if(plotMu){
		if(axes){
			lines(x$x, 0.5*x$meas/x$cdfDist, col=col[3],
					lwd=lwd[3], lty=lty[3])
			axis(4, at=0.5*range(x$meas)/x$cdfDist,
					labels=signif(c(0,x$cdfDist), sigdig))
		} else {
			lines(x$x, x$meas, col=col[3], lwd=lwd[3], lty=lty[3])
		}
		mtext(rlab, 4, line=3)
	}
	if(axes){
		axis(1)
		axis(2)
	}
	if(leg[1] != FALSE){
		legend('topleft', lty=lty, lwd=lwd, col=col, bg='#FFFFFFAA', legend=leg)
	}
}

