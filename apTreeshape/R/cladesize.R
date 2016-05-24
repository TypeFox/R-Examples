"cladesize" <-
function(tree) {
	
	if (class(tree)!='treeshape') {
		stop("invalid arguments")
	}
	
	node=sample(1:nrow(tree$merge),1)
	clade=smaller.clade.spectrum(tree)
	child.number<-clade[node,1]	
	child.number
}

