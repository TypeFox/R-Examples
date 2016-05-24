sister <-
function(phy, node){
	
	if (is.character(node))
		node <- which(phy$tip.label %in% node)
	
	if (node == length(phy$tip.label) + 1) D <- 0
	else {
		x <- phy$edge[,1][phy$edge[,2] == node] # getmrca
		sister <- phy$edge[,2][phy$edge[,1] == x] # desc of mrca
		sister <- sister[sister != node] # eliminate node
		if (sister <= length(phy$tip.label)) D <- sister
		else D <- descendants(phy, sister) # get whole sister clade
	}
	D
}