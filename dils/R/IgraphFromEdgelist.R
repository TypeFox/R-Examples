#' Convert an edgelist to an igraph
#'
#' Given the adjacency matrix for a network return a data.frame listing all possible edges and the weights for each edge.
#' 
#' @param elist data.frame, see 'Details' for formatting assumptions.
#' @param directed logical, If TRUE, the returned igraph is directed.
#' @return igraph, If the edgelist third column has values other than {0, 1} then the weights are stored in E(returned graph)$weight.
#' @export
#' @seealso \code{\link{EdgelistFill}}
#' @references
#' \url{https://github.com/shaptonstahl/}
#' @author Stephen R. Haptonstahl \email{srh@@haptonstahl.org}
#' @details This assumes that \code{elist} is a data.frame with three columns. Each row is an edge in the network. The first column lists the node the edge is coming from, the second column lists the node the edge is going to, and the third column lists the weight of the edge.
#' @examples
#' edgelist <- cbind(expand.grid(letters[1:2], letters[1:2]), runif(4))
#' g <- IgraphFromEdgelist(edgelist)
#' get.edgelist(g)
#' E(g)$weight
#' plot(g, edge.width=5*E(g)$weight, edge.curved=TRUE)
IgraphFromEdgelist <- function(elist, 
                               directed=TRUE) {
  # Guardians
  # elist ones are implemented by AdjacencyFromEdgelist
  stopifnot(is(directed, "logical"))
  
  adj.result <- AdjacencyFromEdgelist(elist)
  g <- graph.adjacency(adj.result$adjacency, 
                       mode=ifelse(directed, "directed", "undirected"),
                       weighted=ifelse(all(adj.result$adjacency %in% c(0,1)), NULL, TRUE))
  V(g)$name <- adj.result$nodelist
  return(g)
}
