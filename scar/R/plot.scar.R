plot.scar <- function(x, whichplot = 1:d, addp = TRUE, ...){
	n=nrow(x$x)
	d=ncol(x$x)
	par(ask=TRUE)
	for(j in whichplot){
		xy=matrix(c((x$x)[,j],(x$componentfit)[,j]),nrow=n, ncol=2)
		sortxy=xy[order(xy[,1]),]
		plot(sortxy[,1],sortxy[,2],type="l",xlab=paste("x",j, sep = ""), ylab="y", main=paste("Component ",j), ...)
		if (addp==TRUE) points(sortxy[,1],sortxy[,2])
	}
	par(ask=FALSE)
}
