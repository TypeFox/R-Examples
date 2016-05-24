diffuse<-function(hexbin, tab, sorted=FALSE){
	n<-nrow(tab)
			
	nclasses<-ncol(tab)	
	slop<-numeric(nclasses)
	allocations<-matrix(nrow=n,ncol=6)	
	xy<-cbind(hexbin@xcm, hexbin@ycm)
	for(i in 1:n){
		alloc<-sainte_lague(tab[i,],6)
		allocations[i,]<-sample(alloc)
                if (sorted) allocations[i,]<-sort(allocations[i,])
		nnearby<-min(6,(n-i))
		slop<-oldslop<-slop+attr(alloc,"error")
		if (nnearby>0){
			nearby<-FNN::knnx.index(xy[-(1:i),,drop=FALSE],xy[i,,drop=FALSE],k=nnearby)
			for(j in nearby) {
				temp<-tab[i+j,]+oldslop/nnearby
				slop<-slop+pmin(0,temp)-oldslop/nnearby
				tab[i+j,]<-pmax(0,temp)
				}
		}
	}	
	t(allocations)
}
