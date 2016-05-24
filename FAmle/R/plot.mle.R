plot.mle <-
function(x,ci=FALSE,alpha=.05,...)
{
	layout(matrix(1:4,2,2))
	Quantile.plot(x,ci,alpha)
	cdf.plot(x)
	Return.plot(x,ci,alpha)
	hist.plot(x)
	title(main=paste('Diagnostic plots : ',x$data.name,' ~ ',x$dist,sep=''),outer=TRUE,line=-2,cex.main=1.5)
	layout(matrix(1,1,1))

}

