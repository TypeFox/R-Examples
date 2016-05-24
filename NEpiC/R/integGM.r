integGM <-
function (G, genes, weights, simplify=T){
	#if (!require(igraph)) {
	#	stop('igraph must be pre-installed!\n')
	#}
	if (is.element("weight", list.vertex.attributes(G))) {
		cat("Warning: previous G node weight replaced!\n")
	}
	names(weights) <- genes
	genes <- intersect(genes,V(G)$name)
	subG <- induced.subgraph(G,genes)
	V(subG)$weight <- weights[V(subG)$name]
	if (simplify) subG <- simplify(subG)
	return(subG)
}

