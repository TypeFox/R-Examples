"clade.members" <-
function(x, phy, tip.labels=FALSE, include.nodes=FALSE){
    
    # NEW2OLD: CONVERTED...
	
	# returns a vector of the tips that descend from an identified node
	if(class(phy) != "phylo") stop("Phylogeny required")
	
	NallNodes <- max(phy$edge)
	Ntips <- max(phy$edge) - phy$Nnode
	
	if(!(x %in% 1:NallNodes)) stop("Node not in range for phylogeny")
	
	# find the children of the node, append them to the vector of nodes (x)
	# and remove the parent, until all the nodes in the vector are tips...
	# now updated to keep track of parents...
	
	intN <- x[x>Ntips]
	descNode <- numeric(length=0)
	
	while(length(intN) > 0){
		
		minIntN <- min(intN)
		childOfMinIntN <- with(phy, edge[,2][which(edge[,1] == minIntN)])
		
		descNode <- c(descNode, minIntN)
		x <- c(x[x != minIntN], childOfMinIntN)
		
		intN <- x[x>Ntips]
	}
	
	RET <- unique(x)
	
	if(tip.labels) {
		RET <- phy$tip.label[x]
	} 
	
	if(include.nodes) {
		RET <- list(tips=RET, nodes=descNode)
	}
	
	return(RET)
	
}

"clade.members.list" <-
function(phy, tips=FALSE, tip.labels=FALSE, include.nodes=FALSE){
    
    # OLD2NEW CONVERTED
    
	# returns a list of vectors showing the tips
	# subtending from each node in the tree
	if(class(phy) != "phylo") stop("Phylogeny required")

	nodes <- 1:max(phy$edge)
	
	if(!tips) nodes <- nodes[nodes > length(nodes) - phy$Nnode]
	
	clade.list <- mapply(clade.members, nodes, MoreArgs=list(phy=phy, tip.labels=tip.labels, include.nodes=include.nodes), SIMPLIFY=FALSE)
	names(clade.list) <- nodes
	
	return(clade.list)
}

## "all.clades" <-
## function(phy, tips=FALSE, tip.labels=FALSE){
##     
##     .Deprecated("clade.members.list")
##     clade.members.list(phy, tips=FALSE, tip.labels=FALSE)
## }
