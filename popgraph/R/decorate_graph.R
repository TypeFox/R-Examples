#' Merge data into popgraph node-wise
#' 
#' This function takes a data frame and a popgraph object and attempts to add
#'  the node data (from the data frame) to the graph as vertex attributes.
#' @param graph An object of type \code{igraph} to populate
#' @param data An object of type \code{data.frame} that has a Node column 
#'  that is the same as the V(graph)$name values.
#' @param stratum The column name of node lables to match up with graph.stratum
#'  (default 'Population') to be matched with V(graph)$name attribute
#' @return A populated igraph object with as much of the metadata in 
#'  the data.frame as possible stitched into the \code{igraph} object
#' @author Rodney J. Dyer <rjdyer@@vcu.edu>
#' @import igraph
#' @export
decorate_graph <- function( graph, data, stratum="Population" ) {
  
  if( missing(graph) | missing(data) ) 
    stop("You need to pass both a graph AND data to the decorate_popgraph function...")
  
  # decorating the node set
  if( is(data,"data.frame")){
    
    if( !(stratum %in% names(data))){
      if( "Stratum" %in% names(data))
        stratum <- "Stratum"
      else
        stop("You must indicate which vertex attribute represents node names.")
    }
    
    
    node.labels <- V(graph)$name 
    num.nodes <- length( node.labels )
    
    data <- data[ (data[[stratum]] %in% node.labels) , ]
    
    #add null data
    cats <- setdiff( names(data), stratum )
    
    # go through the categories and fill them in with null 
    for( cat in cats ){
      vec <- rep( NA, num.nodes )
      for( i in 1:num.nodes ){
        node.name <- node.labels[i]
        
        # only use those 
        if( node.name %in% data[[stratum]] ) {
          val <- data[ data[[stratum]]==node.labels[i],cat ]
          vec[i] <- val        
        }
      }
      graph <- set.vertex.attribute(graph, cat, value=vec )    
    }  
  }
  
  # deocorating the edge set
  else if( is(data, "matrix") ) {
    
    prob <- rep(0,length(E(graph)))
    for( i in 1:length(E(graph))){
      e <- get.edge(graph=graph,id=i)
      prob[i] <- data[e[1],e[2]]
    }
    graph <- set.edge.attribute(graph,"probability",index = E(graph),value = prob)
  }
  
  
  
  return( graph )    

}