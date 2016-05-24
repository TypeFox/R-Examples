distSumConnectMatrix <- function(g, dist=NULL) {
  if (class(g)[1] != "graphNEL")
    stop("'g' has to be a 'graphNEL' object")
  stopifnot(.validateGraph(g))
  
  if (is.null(dist))
    dist <- distanceMatrix(g)

  distdeg <- rowSums(dist)

  e <- edges(g)
  m <- matrix(0, nrow=numNodes(g), ncol=numNodes(g),
              dimnames=list(nodes(g), nodes(g)))
  for (i in names(e))
    for (j in e[[i]])
      m[[i, j]] <- 1 / sqrt(distdeg[[i]] * distdeg[[j]])

  m
}
