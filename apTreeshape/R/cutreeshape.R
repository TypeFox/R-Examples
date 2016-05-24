"cutreeshape" <-
function(tree, node, type) {


bottom.cut <- function (tree, node) {
	clade=smaller.clade.spectrum(tree)
	tip=clade[nrow(tree$merge)-node+1,1]
	
	merge=matrix(NA, tip-1, 2)
	current.node=tip-2
	current.line=tip-1
	
	lines=node
	
	while(identical(lines, numeric(0)) == FALSE){
		merge[current.line,]=tree$merge[lines[1],]
		if (merge[current.line,1]>0) {
			lines=c(lines, merge[current.line,1])
			merge[current.line,1]=current.node
			current.node=current.node-1
		}
		if (merge[current.line,2]>0) {
			lines=c(lines, merge[current.line,2])
			merge[current.line,2]=current.node
			current.node=current.node-1
		}
		lines=lines[-1]
		current.line=current.line-1
	}
	
	names=c()
	current.tip=-1
	for (i in 1:length(merge)) {
		if (merge[i]<0) {
			names=c(names, merge[i])
			merge[i]=current.tip
			current.tip=current.tip-1
		}
	}
	names=-names
	names=tree$names[names]
	res=treeshape(merge, names)
	res 
}


top.cut <- function(tree, node) {
#We keep the top part of the tree
	merge=matrix(nrow=nrow(tree$merge)-node, ncol=2)
	for (i in 1:nrow(merge)){
		merge[i,]=tree$merge[node+i,]
	}
	currenttip=-1
	for (i in 1:length(merge)){
		if (merge[i]<=node) {
			merge[i]=currenttip
			currenttip=currenttip-1
		} else {
			merge[i]=merge[i]-node
		}
	}
	res=treeshape(merge)
	res
}



	if (class(tree)!='treeshape') {
		stop("invalid arguments1")
	}
	if (node!=floor(node) | node < 0){
		stop("node must be a positive integer")
	}
	
	if (type=="top"){
		return(top.cut(tree, node))
	}
	if (type=="bottom"){
		return(bottom.cut(tree, node))
	}
	stop("wrong values for parameter type")
	

}

