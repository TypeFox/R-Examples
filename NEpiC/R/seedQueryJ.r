seedQueryJ <-
function (G, seed, search_r = 2, r = 0.1){
	#if (!require(igraph)) {
	#	stop('igraph must be pre-installed!\n')
	#}
	net <- G
	d <- search_r
	if (!is.element("name", list.vertex.attributes(net))) {
		stop("Graph must have 'name' attribute")
	}
	if (!is.element("weight", list.vertex.attributes(net))) {
		stop("Graph must have 'weight' attribute")
	}
	subG <- induced.subgraph(net, seed)
	if (!is.connected(subG)) {   		#the seed must be connected
		stop("Input seeds are disjoint")
	}
	in.nodes <- V(subG)$name
	while (TRUE) {
		subx <- V(subG)$name
		for (rad in 1:d) {
			subsum <- sum(V(subG)$weight)/sqrt(length(subx)) #Peilin
			
			tmp.neigh <- unlist(neighborhood(net, order = rad, nodes = V(subG)$name)) 
			pot.nodes <- V(net)[tmp.neigh]$name
			pot.nodes <- setdiff(pot.nodes, in.nodes)
			if (length(pot.nodes) == 0) break
			sub.weg <- V(net)[pot.nodes]$weight
			best.nodes <- pot.nodes[which(sub.weg == max(sub.weg))]
			
			subsum.u <- (sum(V(subG)$weight) + V(net)[best.nodes[1]]$weight)/sqrt(length(subx)+1)
			
			if (subsum.u > subsum * (1 + r)) {
				tmp <- unlist(lapply(best.nodes, function(x) node2treePath(net,V(subG)$name, x)))
				in.nodes <- c(tmp, V(subG)$name)
				subG <- induced.subgraph(net, in.nodes)
				break
			}
		}
		if (length(subx) == vcount(subG)) break
	}
	return(subG)
}

