node2treePath <-
function (G, Tnodes, node){
	#if (!require(igraph)) {
#		stop('igraph must be pre-installed!\n')
#	}
	tmp.path <- get.all.shortest.paths(G, node, Tnodes)$res
	tmp.l <- unlist(lapply(tmp.path, length))
	index <- which(tmp.l == min(tmp.l))
	
	tmp.path = tmp.path[index]
	tmp.sum <- unlist(lapply(tmp.path, function(x)return(sum(V(G)[x]$weight))))
	index <- which(tmp.sum == max(tmp.sum))
	
	selected.path = tmp.path[index]
	collect <- unlist(lapply(selected.path, function(x)return(V(G)[x]$name)))
	
	return(collect)
}

