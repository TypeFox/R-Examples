tree.na.Count <-
function(tree,thr=6,srt="drop",min.val=0.01,mod.id=c(1,0,0,0)){
	N=length(tree$edge.l)
	subs=subtrees2.6(tree)
	X <- 0
	for(i in 1:(N/2)){
	sZ<-length(subs[[i]]$edge.length)
		if(sZ<=thr){X <- X+1}
		if(sZ>=(N-thr)){X <- X+1}
		}
	X
	}

