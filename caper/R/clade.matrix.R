"clade.matrix" <-
function(phy){

    # OLD2NEW: CONVERTED
    
	# returns the phylogeny as a table showing clade
	# membership below each node. 

	# check the object is a phylogeny
	if(class(phy) != "phylo") stop("Phylogeny required")
	
	# get the number of tips and nodes
	nb.tips <- max(phy$edge) - phy$Nnode
	nb.nodes <- phy$Nnode
	nb.all <- nb.tips + nb.nodes
	
	# set-up the clade matrix
	mat.names <- list(edges= c(1:nb.all), tips=1:nb.tips)
	clade.mat <- matrix(0, nrow=nb.all, ncol=nb.tips, dimnames=mat.names)
	
	# the diagonal from [1,1] are the tips
	diag(clade.mat) <- 1
	
	# now deal with internals
	node.members <- clade.members.list(phy)

	node.id <- names(node.members)
	for(rows in seq(along=node.id)){
		clade.mat[node.id[rows], node.members[rows][[1]]] <- 1}

	RET <- list(clade.matrix=clade.mat, tip.label=phy$tip.label, edge=phy$edge)
	class(RET) <- "clade.matrix"

	# if they exist, get edge lengths into correct order, inserting root edge
	if(! is.null(phy$edge.length)){
		if(is.null(phy$root.edge)){
			edge.len <- c(phy$edge.length,0)
			} else {
			edge.len <- c(phy$edge.length,phy$root.edge)}
	
		names(edge.len) <- c(phy$edge[,2], nb.tips + 1)
		edge.len <- edge.len[as.character(mat.names$edges)]
		RET$edge.length <- edge.len
	}
	
	return(RET)
}

