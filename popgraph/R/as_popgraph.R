#' Converts an object of type matrix or igraph into \code{popgraph}
#' 
#' This is a simple conversion routine for \code{matrix} or \code{igraph} 
#'  objects into \code{popgraph} objects
#' @param graph An object of type \code{matrix} or \code{igraph}
#' @return An object of type \code{popgraph}
#' @export
#' @author Rodney J. Dyer <rjdyer@@vcu.edu>
as.popgraph <- function(graph) {
  ret <- NULL
  if( is(graph,"matrix")) {  
    
    ret <- igraph::graph.adjacency( graph, mode="undirected",weighted=TRUE) 
    if( is.null(colnames( graph )) )
      V(ret)$name <- as.character(paste("node",seq(1,ncol(graph)), sep="-"))
    else
      V(ret)$name <- colnames( graph )
  }
  
  if( is(graph,"igraph")) {
    ret <- graph
    if( length(igraph::V(ret)$name) == 0)
      igraph::V(ret)$name <- paste( "Node",seq(1,length(V(ret))),sep="-")
  }
    
  
  if( is.null(ret) )
    stop("Must pass either an igraph or matrix item to this function.")
  
  class(ret) <- c("popgraph","igraph")
  return(ret)
    
}