net.edges <-
function(theta)
{
	library(igraph)
	adj = make.adj.matrix(theta,separate=TRUE)
	K = length(theta)

	edges = list()
	for(k in 1:K)
	{
		diag(adj[[k]])=0
		gadj = graph.adjacency(adj[[k]],mode="upper")
		edges[[k]] = E(gadj)
	}
	return(edges)
}

