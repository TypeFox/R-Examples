informationLayerIndex <- function(g, dist=NULL, layer=NULL) {
  if (class(g)[1] != "graphNEL")
    stop("'g' must be a 'graphNEL' object")
  if (is.null(dist))
    dist <- distanceMatrix(g)
  if (is.null(layer))
    layer <- layerMatrix(g, dist=dist)

  nn <- numNodes(g)

  lis <- sapply(nodes(g), function(i) {
    sum(sapply(1:max(dist[i,]), function(j) {
      p <- layer[[i, j]] / nn
      p * log2(p)
    }))
  })

  -sum(lis)
}
