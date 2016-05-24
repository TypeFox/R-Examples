"all.equal.treeshape" <-
function(target, current,names=FALSE, height=FALSE, ...){
	tree1 <- target
	tree2 <- current

	clade1 <- smaller.clade.spectrum(tree1)
	clade2 <- smaller.clade.spectrum(tree2)
	
aux <- function(node1, node2) {
	if ((node1<0) & (node2<0)) {
		return((tree1$names[-node1]==tree2$names[-node2]) | (!names))
	}
	if ((node1*node2)<0) {return(FALSE)}
	
	if (identical(clade1[node.number-node1+1,],clade2[node.number-node2+1,])==FALSE) {
		return(FALSE)
	}
	
	res1=aux(tree1$merge[node1,1], tree2$merge[node2,1])
	res2=aux(tree1$merge[node1,2], tree2$merge[node2,2])
	res3=aux(tree1$merge[node1,2], tree2$merge[node2,1])
	res4=aux(tree1$merge[node1,1], tree2$merge[node2,2])
	
	if (height) {
		aux1=res1&res2&(tree1$merge[node1,1]==tree2$merge[node2,1])&(tree1$merge[node1,2]==tree2$merge[node2,2])
		aux2=res3&res4&(tree1$merge[node1,1]==tree2$merge[node2,2])&(tree1$merge[node1,2]==tree2$merge[node2,1])
	} else {
		aux1=res1&res2
		aux2=res3&res4
	}
	return(aux1 | aux2)
}

	if (nrow(tree1$merge)==nrow(tree2$merge)) {
		node.number=nrow(tree1$merge)
	} else {
		return(FALSE)
	}

	res=aux(nrow(tree1$merge), nrow(tree2$merge))
	res
}

