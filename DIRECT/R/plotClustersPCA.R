plotClustersPCA <-
function (item.names, data.summary,PCA.label.adj=-0.01,...)
{
	nclust = data.summary$nclust
	clust.labels = data.summary$top.clust.labels
	clust.ind = data.summary$top2allocations[,1]
	post.probs = data.summary$post.alloc.probs
	
	clust.PCA = princomp (post.probs, scores=TRUE)
	print (clust.PCA)
	
	col.tmp = rainbow(nclust)
#	par (mfrow=c(1,1), mar=c(4.5, 4.5, 0.5, 0.5))
	plot (clust.PCA$scores[,1:2], xlab="PC1", ylab="PC2", pch=16, cex=0.4)
	for (i in 1:nclust)
		text (clust.PCA$scores[which (clust.ind==clust.labels[i]),1], clust.PCA$scores[which (clust.ind==clust.labels[i]),2]+PCA.label.adj, labels=item.names[which (clust.ind==clust.labels[i])], cex=0.5, col=col.tmp[i])
	abline (h=0, lty=2)
	abline (v=0, lty=2)
	
}

