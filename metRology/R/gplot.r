#Provides a simple grouped plot of a matrix using mode "h"
#A plot of values in x grouped by row 
#Called by plot.mandel

gplot <- function(x, main=NULL, xlab=NULL, ylab=deparse(substitute(x)) ,
				ylim=NULL, las=1, axes=TRUE, cex.axis=1, frame.plot = axes,
				lwd=1, lty=1,col=par("col"), 
				separators=TRUE, col.sep="lightgrey", lwd.sep=1, lty.sep=1,
				zero.line=TRUE, lwd.zero=1, col.zero=1, lty.zero=1, ...) {
		
		if(length(dim(x)) != 2) stop("gplot requires an object with two dimensions")
		if(is.null(row.names(x))) stop("gplot requires an object with valid row names")
		
		pars<-list(...)
		
		ni<-ncol(x)
			#Number of items to plot per group
		
		ng <- nrow(x)
			#Number of groups
				
		stkx<-stack(as.data.frame(t(x)))
		
		levels(stkx$ind) <- row.names(x)
			#Returns level order to original row name order
			
		mids <- 1:ng
		
		sep <- min(0.3, 0.6/max(1, ni-1))
		
		offsets <- cumsum(rep(sep, ni))
		
		offsets <- offsets - mean(offsets)
		
		locations <- rep(mids, each=ni) + rep(offsets, ng)
		
		if(is.null(ylim)) ylim <- range(pretty(c(0, na.omit(stkx$values))))
		
		plot(locations, stkx$values, type="n", 
			ylim=ylim, main=main, xlab=xlab, ylab=ylab, axes=FALSE, 
			col=col, ...)
		if(separators) abline(v=c(0, mids)+0.5, col=col.sep, lty=lty.sep, lwd=lwd.sep)
		
		segments(locations, rep(0, length(locations)), locations, stkx$values, 
			lwd=lwd, col=col, lty=lty)	

		if(zero.line) abline(h=0, col=col.zero, lwd=lwd.zero, lty=lty.zero)

		
		if(axes) {
			axis(2, las=las, cex.axis=cex.axis)
			axis(1, at=mids, labels=levels(stkx$ind), las=las, cex.axis=cex.axis)
		}
		
		if(frame.plot) box()
		
		
		return(invisible(mids))
}

