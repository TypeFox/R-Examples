plotClustersMean <-
function (data, data.summary, SKIP, nTime=length(times), times=1:nTime,...)
{
	if (nrow(data) != data.summary$nitem)
		stop ("data and summary statistics do not match")
	
	nItem = nrow (data)
	nRep = (ncol (data) - SKIP) / nTime
	
	ts = array (0, dim = c(nItem, nTime, nRep))
	for (r in 1:nRep)
	{
		ts[,,r] = as.matrix (data[,SKIP + (0:(nTime-1))*nRep + r])
	}
	
	ts.mean = apply (ts, c(1,2), mean)
	
	nclust = data.summary$nclust
	clust.labels = data.summary$top.clust.labels
	post.clust.pars.mean = data.summary$post.clust.pars.mean
	clust.ind = data.summary$top2allocations[,1]
	
	col.tmp = rainbow(nclust)
	plot.nrow = round (sqrt (nclust))
	plot.ncol = round (nclust / plot.nrow)
	if (plot.nrow*plot.ncol < nclust)	plot.ncol=plot.ncol+1
	par (mfrow=c(plot.nrow, plot.ncol), ...)
	for (i in 1:nclust)
	{
		tmp = which (clust.ind==clust.labels[i])
		matplot (times, t(ts.mean[tmp,]), xlab="", ylab="", main=paste("Cluster ", i, " (", length(tmp), " genes)", sep=""), type="l", lty=1, col=1, lwd=3)
		lines (times, post.clust.pars.mean[i,1:nTime], lty=1, lwd=3, col=col.tmp[i])
		abline (h=0, lty=2, col="brown", lwd=2)
	}
	
}

