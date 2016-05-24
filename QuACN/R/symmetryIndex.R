symmetryIndex <- function(g, dist=NULL, deg=NULL) {
  if (class(g)[1] != "graphNEL")
    stop("'g' has to be a 'graphNEL' object")
  stopifnot(.validateGraph(g))
  
  if (is.null(dist))
    dist <- distanceMatrix(g)
  if (is.null(deg))
    deg <- graph::degree(g)

  ig <- .G2IG(g)
  Ai <- .cardNi(g, dist, deg)
  aut <- as.numeric(igraph::graph.automorphisms(ig)[["group_size"]])

  1/numNodes(g) * sum(Ai * log2(Ai)) + log2(aut)
}

#
# .orbitSizes <- function(ig, dist) {
#   vertices <- c(1:vcount(ig))
#   Ai <- c()
#
#   while (length(vertices) != 0) {
#     focus <- vertices[1]
#     i <- 2
#     while (i <= length(vertices)) {
#       current <- vertices[i]
#       if (max(dist[focus,]) == max(dist[current,])) {
#         sapply(1:max(dist[focus,]), function(j) {
#         })
#       }
#       i <- i + 1
#     }
#     vertices <- vertices[-1]
#   }
#
#   Ai
# }
