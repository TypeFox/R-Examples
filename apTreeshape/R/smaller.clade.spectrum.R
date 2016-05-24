"smaller.clade.spectrum" <-
function(tree) {
	if (identical(tree,NULL)) {
		stop("invalid tree","\n")
	}
	merge=tree$merge
	
	spec <- matrix(0,nrow=nrow(merge),ncol=2)
	spec[nrow(merge),]=c(2,1)
	
	for (node in 2:nrow(merge)) {
		if (merge[node,1]<0) {lc=1}
		else {lc=spec[nrow(merge)-merge[node,1]+1, 1]}
		
		if (merge[node,2]<0) {rc=1}
		else {rc=spec[nrow(merge)-merge[node,2]+1, 1]}
		
		spec[nrow(merge)-node+1,]=c(rc+lc, min(rc, lc))
		
	}
	spec
	
}

