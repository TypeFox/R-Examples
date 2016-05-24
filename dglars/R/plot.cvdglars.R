plot.cvdglars <- function(x,...){
	dev_m <- x$dev_m
	dev_v <- x$dev_v
	nfold <- x$control$nfold
	k <- qt(0.975,nfold-1)
	dev_up <- dev_m + k * sqrt(dev_v/nfold)
	dev_low <- dev_m - k * sqrt(dev_v/nfold)
	ng <- x$control$ng
	g <- seq(x$g_max,x$g0,length=ng)
	g_hat <- g[which.min(dev_m)]
	df <- sum(abs(x$beta)>0)
	plot(g,dev_m,xlab=expression(gamma),ylab="Deviance",ylim=c(min(dev_low),max(dev_up)),pch=20,type="n",main="Cross-Validation Deviance",...)
	segments(x0=g,y0=dev_low,y1=dev_up,col=8,lty=2)
	points(g,dev_low,,pch="-")
	points(g,dev_m,pch=20)
	points(g,dev_up,,pch="-")
	abline(v=g_hat,col=2,lty=2,lwd=2)
	axis(3,at=g_hat,labels=paste("df = ",df,sep=""),padj=1)
}
