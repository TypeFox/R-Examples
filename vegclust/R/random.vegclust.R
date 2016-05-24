random.vegclust<-function(x, cmin=2, cmax=20, nstart=10, verbose=TRUE,...) {
	nvc = cmax-cmin +1
	vc = vector("list",nvc)
	i=1
	for(c in cmin:cmax) {
	   if(verbose) cat(paste("PROCESSING",c,"MOBILE CLUSTERS\n"))
	   vc[[i]] = vegclust(x,mobileCenters=c,nstart=nstart,...)
		i = i+1
	}
	class(vc)<-"mvegclust"
	return(vc)
}