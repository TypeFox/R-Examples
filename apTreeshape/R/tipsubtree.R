"tipsubtree" <-
function (tree, tips, numeric=FALSE) {
if (identical(tree,NULL)) {
		stop("invalid tree","\n")
		
}		
	tip <- c()
	if (!numeric){
		for (i in tips){ tip=c(tip,which(tree$names == i)) }
	}
	else {
		tip=tips
	}
	tip=sort(tip)
	if (length(tip)<2){
		stop("at least two tips")
	}
	
	na=(1:length(tree$names))[-tip]
	na=c(0, na)
	na=-na
	
	for (i in 1:nrow(tree$merge)) {
		tmp=0
		if((sum(tree$merge[i,1]==na))==1){
			tree$merge[i,1]=0
			tmp=tmp+1
		}
		if((sum(tree$merge[i,2]==na))==1){
			tree$merge[i,2]=0
			tmp=tmp+2
		}
		
		if (tmp==1){
			tree$merge[which(tree$merge==i)]=tree$merge[i,2]
			tree$merge[i,2]=0
		}
		if (tmp==2){
			tree$merge[which(tree$merge==i)]=tree$merge[i,1]
			tree$merge[i,1]=0
		}
		if (tmp==3){
			tree$merge[which(tree$merge==i)]=0
		}
	
	}
	merge=matrix(NA, nrow=(length(tip)-1), ncol=2)
	
	nodes=list()
	current.node=1
	current.tip=-1
	names=c()
	
	for (i in 1:nrow(tree$merge)) {
		if (tree$merge[i,1]*tree$merge[i,2]==0){
#			na.node=na.node+1
		} else {
			merge[current.node, ]=tree$merge[i,]
			nodes[[i]]=current.node
			
			if (merge[current.node,1]>0){
				merge[current.node,1]=nodes[[merge[current.node,1]]]
			} else {
				names=c(names, tree$names[-merge[current.node,1]])
				merge[current.node,1]=current.tip
				current.tip=current.tip-1
			}
			if (merge[current.node,2]>0){
				merge[current.node,2]=nodes[[merge[current.node,2]]]
			} else {
				names=c(names, tree$names[-merge[current.node,2]])
				merge[current.node,2]=current.tip
				current.tip=current.tip-1
			}
			current.node=current.node+1
		}
	
	}
	res=treeshape(nodes=merge, names=names)
	
	res
}

