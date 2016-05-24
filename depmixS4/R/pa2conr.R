
# convert conpat vector to rows of constraint matrix

pa2conr <-
function(x) {
	fix=as.logical(x)
	x=replace(x,list=which(x==1),0)
	un=setdiff(unique(x),0)
	y=matrix(0,0,length(x))
	for(i in un) {
		z=which(x==i)
		for(j in 2:length(z)) {
			k=rep(0,length(x))
			k[z[1]]=1
			k[z[j]]=-1
			y=rbind(y,k)
		}
	}
	pa = list(free=fix,conr=y)
	return(pa)
}

