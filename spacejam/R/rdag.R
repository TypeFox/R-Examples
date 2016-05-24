rdag <- function (p, nedges) 
{
    stopifnot(nedges >= 0, nedges <= choose(p,2), as.integer(nedges) == nedges)
	M <- matrix(0,p,p)
    M[upper.tri(M)][sample(choose(p,2),nedges)] <- 1
    graph.adjacency(adjmatrix = M)
}