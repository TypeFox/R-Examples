frEZ <- function(y, d = 2){ # from library(igraph)
	x <- layout.fruchterman.reingold(graph.adjacency(y), dim = d)
	x / apply(x , 2, sd)
}