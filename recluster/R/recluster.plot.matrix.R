recluster.plot.matrix<-function(mat)
{
	dist<-as.matrix(mat)
	plot(c(1, nrow(dist)), c(1, nrow(dist)), type = "n")
	for (i in 1:nrow(dist)){
		for (n in 1:ncol(dist)){
		cols <- grey((dist[i,n] - min(dist))/(max(dist)-min(dist)))
		rect(nrow(dist)-i, n, nrow(dist)-1-i, 1+n, col=cols, 		border = NA , lwd=0)
		}
	}
}
