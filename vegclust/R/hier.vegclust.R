hier.vegclust<-function(x, hclust, cmin=2,cmax=20, verbose=TRUE,...) {
	vc = vector("list",length=cmax-cmin +1)
	xs = x[row.names(x)%in% hclust$labels,] #Select those objects used for hierarchical clustering
	i=1
	for(c in cmin:cmax) {
	   if(verbose) cat(paste("PROCESSING",c,"MOBILE CLUSTERS\n"))
	   cent = as.vegclust(xs,cutree(hclust,k=c))$mobileCenters
		vc[[i]] = vegclust(x,mobileCenters=cent,...)
		i = i+1
	}
	class(vc)<-"mvegclust"
	return(vc)
}