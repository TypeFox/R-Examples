globalQueryJ <-
function (G, search_r = 2, r = 0.2, min.size=5){
#	if (!require(igraph)) {
#		stop('igraph must be pre-installed!\n')
#	}	
	sublist = list()
	for (node in V(G)$name) {
		ng <- seedQueryJ(G, node, search_r, r)
		if (vcount(ng) >= min.size) sublist[[node]] <- ng       #minimum size filtering
	}
	return(sublist)
}

