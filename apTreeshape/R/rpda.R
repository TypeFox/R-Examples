"rpda" <-
function(tip.number) {
	
reorg.iterative <- function (merge, root) {

	node=list()
	current.node=nrow(merge)
	node[[root]]=current.node
	current.node=current.node-1

	while(identical(root, numeric(0)) == FALSE){
		if (merge[root[1],1]>0) {
			root=c(root, merge[root[1],1])
			node[[merge[root[1],1]]]=current.node
			current.node=current.node-1
		}
		if (merge[root[1],2]>0) {
			root=c(root, merge[root[1],2])
			node[[merge[root[1],2]]]=current.node
			current.node=current.node-1
		}
		root=root[-1]
	}
	
	final.merge=merge
	
	for (i in 1:nrow(merge)) {
		n=node[[i]]
		final.merge[n,]=merge[i,]
		if (final.merge[n,1]>0) {
			final.merge[n,1]=node[[final.merge[n,1]]]
		}
		if (final.merge[n,2]>0) {
			final.merge[n,2]=node[[final.merge[n,2]]]
		}
		
	}
	
	final.merge

}

	
			
	if (tip.number<2 | tip.number!=floor(tip.number)) {
		stop("tip.number must be an integer greater than 2")
	}

	if (tip.number==2) {
		tmp=c(-1,-2)
		res=list(merge=matrix(1, ncol=2, nrow=1))
		res$merge[1,]=c(-1,-2)
		class(res)<-'treeshape'
		res
	
	}
	else{
		merge<-matrix(1, ncol=2, nrow=(tip.number-1))
		merge[1,]=c(-1,-2)
		
		root=1
		
		for (tip in 3:tip.number) {
			
			n<-floor(runif(n=1, min=1, max=2*tip-2))
			if (n==1) {
				tmp=c(root, -tip)
				root=tip-1
				merge[tip-1,]=tmp
			}
			else {
				div=floor(n/2)
				mod=n-2*floor(n/2)
				
				tmp=merge[div,mod+1]
				merge[div,mod+1]=tip-1
				merge[tip-1,]=c(tmp, -tip)
			}
		}
		
		final.merge=reorg.iterative(merge, root)
		res=treeshape(final.merge)
		res
	}
	
}

