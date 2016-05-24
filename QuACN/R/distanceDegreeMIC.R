distanceDegreeMIC <- function(g, dist=NULL) {
  if (class(g)[1] != "graphNEL")
    stop("'g' must be a 'graphNEL' object")
  if (is.null(dist))
    dist <- distanceMatrix(g)

  distdeg <- rowSums(dist)
  W <- wiener(g, dist=dist)
  pis <- distdeg / 2 / W

  -sum(pis * log2(pis))
}
