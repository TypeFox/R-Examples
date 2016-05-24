plot.scair <- function(x, whichplot = 1:m, addp = TRUE, ...){
	n=nrow(x$x)
	m=length(x$shape)
	par(ask=TRUE)
	for(j in whichplot){
		xy=matrix(c((x$x %*% x$index)[,j],(x$componentfit)[,j]),nrow=n, ncol=2)
		sortxy=xy[order(xy[,1]),]
		plot(sortxy[,1],sortxy[,2],type="l",xlab=paste("A^T x",j, sep = ""), ylab="y", main=paste("Index ",j), ...)
		if (addp==TRUE) points(sortxy[,1],sortxy[,2])
	}
	par(ask=FALSE)
}
