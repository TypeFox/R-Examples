#' Create a "Network" Object from an igraph Object
#'
#' This function will convert an igraph object to an object compatible with
#' snowboot functions.
#'
#' @param in_graph An igraph object. To create igraph objects from field data,
#'      see \code{\link[igraph]{graph_from_edgelist}},
#'      \code{\link[igraph]{graph_from_data_frame}},
#'      \code{\link[igraph]{graph_from_adjacency_matrix}}, or
#'      \code{\link[igraph]{read_graph}}.
#' @references \url{http://igraph.org/}
#' @return A list that contain elements:
#'    \item{edges}{The edgelist of the network. A two column
#'      \code{matrix} where each row is an edge.}
#'    \item{degree}{The degree sequence of the network, which is
#'      an \code{integer} vector of length n.}
#'    \item{n}{The network order.}
#' @export
#' @examples
#' hex_ring <- igraph::make_ring(6, directed = FALSE, mutual = FALSE, circular = TRUE)
#' net <- igraph_to_network(hex_ring)
igraph_to_network <- function(in_graph){
      if(igraph::is.directed(in_graph))
            stop('Only undirected graphs are supported at the moment.
                  Please import undirected igraph object, e.g. output of
                  graph_from_edgelist(el,directed = F)')
      edges <- igraph::as_edgelist(in_graph)
      edges <- order.edges(edges, ord.col = TRUE)
      degree <- igraph::degree(in_graph)
      #degree.left <- rep(0,igraph::gorder(in_graph))
      n <- igraph::gorder(in_graph)
      list(edges = edges, degree = degree, n = n)
}
