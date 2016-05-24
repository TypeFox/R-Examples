"as.treeshape.treebalance" <-
function(x, ...) {
	tree=x
	
	height=nrow(tree$merge)
	merge=matrix(NA, height, 2)
	
	current.tip=-1
	
	for (node in 1:height) {
		new.node=height+1-node
		
		if (tree$merge[node,2]==1) {
			res1=current.tip
			current.tip=current.tip-1
		} else {
			res1=new.node-1
		}
		
		if ((tree$merge[node,1]-tree$merge[node,2])==1) {
			res2=current.tip
			current.tip=current.tip-1
		} else {
			res2=new.node-tree$merge[node,2]
		}
		merge[new.node,]=c(res1,res2)
	}
	res=treeshape(merge)
}

