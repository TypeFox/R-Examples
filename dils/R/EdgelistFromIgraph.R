#' Convert an igraph to filled edgelist
#'
#' Given an igraph object for a network return a data.frame listing all possible edges and the weights for each edge.
#' 
#' @param g igraph, from \link{igraph} package.
#' @param useWeight logical, Should E(g)$weight be used as the weights for the edges?
#' @return data.frame, full list of all possible edges with weights for each in third column.
#' @export
#' @seealso \code{\link{EdgelistFromAdjacency}}
#' @references
#' \url{https://github.com/shaptonstahl/}
#' @author Stephen R. Haptonstahl \email{srh@@haptonstahl.org}
#' @details This function is preferred to the igraph function \code{get.edgelist} because \code{get.edgelist} only returns rows for edges that have non-zero weight and does not return weights, if present.
#' @examples
#' g <- erdos.renyi.game(10, 2/10)
#' EdgelistFromIgraph(g)
#' 
#' V(g)$name <- letters[1:vcount(g)]
#' EdgelistFromIgraph(g)
#' 
#' E(g)$weight <- runif(ecount(g))
#' EdgelistFromIgraph(g, useWeight=TRUE)
EdgelistFromIgraph <- function(g, 
                               useWeight=FALSE) {
  if(useWeight) {
    A <- get.adjacency(g, attr="weight")
  } else {
    A <- get.adjacency(g)
  }
  if( is.null(V(g)$name) ) {
    node.labels <- paste("node", 1:vcount(g), sep="")
  } else {
    node.labels <- V(g)$name
  }
  out <- EdgelistFromAdjacency(A, nodelist=node.labels)
  return(out)
}
