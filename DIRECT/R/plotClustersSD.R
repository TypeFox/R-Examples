plotClustersSD <- function (data.summary, nTime, ...) {
	nclust <- data.summary$nclust
	post.clust.pars.med <- data.summary$post.clust.pars.median

	col.tmp <- rainbow(nclust)
	plot (1:nclust, post.clust.pars.med[,nTime+1], ylim=c(0,2*max(post.clust.pars.med[, nTime+1:3])), xaxt="n", xlab="cluster", ylab="posterior median of standard deviation", main="", type="p", col=col.tmp, pch=16, cex=1.5, cex.axis = 1.6, cex.lab = 1.6, cex.main=1.6)
	axis (side=1, at=1:nclust, labels=1:nclust, cex.axis=1.6, cex.lab=1.6)
	legend ("topright", legend = c ("within cluster", "across time points", "between replicates"), pch=c(1,2,5), cex=1.5)
	abline (v=1:nclust, lty=2, col=1)
	points (1:nclust, post.clust.pars.med[,nTime+2], pch=17, col=col.tmp, cex=1.5)
	points (1:nclust, post.clust.pars.med[,nTime+3], pch=18, col=col.tmp, cex=1.5)	
}

