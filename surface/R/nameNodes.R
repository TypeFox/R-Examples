nameNodes <-
function(tree){
	if(!is.null(tree$node.label)){
		if(length(unique(tree$node.label))==length(tree$node.label)){
			return(tree)
			}}
	tree$node.label<-paste("zzz",1:tree$Nnode,sep="")
	tree
	}
