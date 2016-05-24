layerMatrix <- function(g, dist=NULL) {
  if (class(g)[1] != "graphNEL")
    stop("'g' must be a 'graphNEL' object")
  stopifnot(.validateGraph(g))
  
  if (is.null(dist))
    dist <- distanceMatrix(g)

  diam <- max(dist)
  t(sapply(nodes(g), function(v) {
    tabulate(dist[v,], nbins=diam)
  }))
}
