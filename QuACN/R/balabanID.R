balabanID <- function(g, dist=NULL) {
  if (class(g)[1] != "graphNEL")
    stop("'g' has to be a 'graphNEL' object")
  stopifnot(.validateGraph(g))

  if (is.null(dist))
    dist <- distanceMatrix(g)

  distdeg <- rowSums(dist)

  n <- nodes(g)
  e <- edges(g)
  .weightedPathSum(e, NULL, function(i, from, to) (1 / sqrt(distdeg[[from]] * distdeg[[to]])))
}
