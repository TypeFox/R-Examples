#' Returns congruence topology
#' 
#' Takes two graphs and returns topology that is 
#'  the intersection of the edge sets.
#' @param graph1 An object of type \code{popgraph}
#' @param graph2 An object of type \code{popgraph}
#' @param warn.nonoverlap A flag indicating that a warning should be thrown
#'  if the node sets are not equal (default = TRUE )
#' @return An object of type \code{popgraph} where the node and edge 
#'  sets are the intersection of the two.
#' @author Rodney J. Dyer <rjdyer@@vcu.edu>
#' @export
congruence_topology <- function( graph1, graph2, warn.nonoverlap=TRUE ) {
  
  if( !inherits(graph1, "popgraph") | !inherits(graph2, "igraph") )
    stop("congruence.topology() requires that you pass an igraph or popgraph object.")
  
  if( warn.nonoverlap & length(setdiff( V(graph1)$name, V(graph2)$name )))
    warning("These two topologies have non-overlapping node sets!  Careful on interpretation.")
  
  nodes <- V(graph1)$name
  for( i in 1:length(nodes)){
    if( !(nodes[i] %in% V(graph2)$name))
      nodes[i] <- NA
  }
  cong.nodes <- nodes[ !is.na(nodes) ]
  #cong.nodes <- intersect( V(graph1)$name, V(graph2)$name )
  a <- as.matrix( get.adjacency( induced.subgraph( graph1, cong.nodes ) ) )
  nms <- row.names(a)
  b <- as.matrix( get.adjacency( induced.subgraph( graph2, cong.nodes ) ) )
  b <- b[nms,nms]
  
  #a <- as.matrix( get.adjacency(graph1))
  #b <- as.matrix( get.adjacency(graph2))
  
  cong <- graph.adjacency( a*b, mode="undirected" )
  class(cong) <- c("igraph","popgraph")
  return( cong )

}