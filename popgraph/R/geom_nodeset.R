#' Plotting of a population graph edge set using ggplot neumonic
#' 
#' This function allows you to layer the edgeset of a \code{popgraph}
#'  object
#' @param mapping The aesthetic mapping as an \code{aes()} object.  This aesthetic
#'  must at least have values for x and y but you can also specify color, 
#'  fill, alpha, label, and size.  Color, alpha, and fill may be a single value
#'  for all nodes or vertex attributes for each node.  
#' @param graph The popgraph/igraph object to be plot
#' @param ... Largely ignored.
#' @return A formatted geom_segment object for addition to a ggplot()
#' @author Rodney J. Dyer <rjdyer@@vcu.edu>
#' @export
#' @examples
#' library(igraph)
#' a <- matrix( c(0,1,0,1,1,0,0,1,0,0,0,1,1,1,1,0),nrow=4)
#' rownames(a) <- colnames(a) <- LETTERS[1:4]
#' graph <- as.popgraph(a)
#' igraph::V(graph)$x <- runif(4)
#' igraph::V(graph)$y <- runif(4)
#' require(ggplot2)
#' ggplot() + geom_nodeset( aes(x=x,y=y), graph )
#' igraph::V(graph)$group <- c("A","A","B","B")
#' ggplot() + geom_nodeset( aes(x=x,y=y,color=group), graph, size=4 )
geom_nodeset<- function( mapping=NULL, graph=NULL, ... ) {
  
  # catch errors with missing 
  if( is.null(mapping))
    stop("You need at least aes(x,y) for aesthetic mapping in this function.")
  
  if( is.null(graph))
    stop("You cannot plot a graph without a graph...")
  
  if( is.null(mapping$x) | is.null(mapping$y))
    stop("To plot a graph, you need coordinates and they must be attributes of the vertices in the graph.")

  # grab mapping labels not in the vertex attributes
  vertex.attr <- list.vertex.attributes(graph)
  mappingNames <- names(mapping)
  for( name in mappingNames) {
    key <- as.character(mapping[[name]])
    if( !(key %in% vertex.attr))
      stop(paste("Aesthetic mapping variable ",key," was not found in the vertex attributes of this graph",sep=""))
  }
  
  x <- y <- color <- colour <- size <- NULL
  
  df <- data.frame( x=get.vertex.attribute(graph,mapping$x), 
                    y=get.vertex.attribute(graph,mapping$y),
                    colour="black")
  
  # size and color
  if( !is.null(mapping$colour) & !is.null(mapping$size)){
    df$color <- get.vertex.attribute(graph,mapping$colour)
    df$size <- get.vertex.attribute(graph,mapping$size)
    ret <- geom_point( aes(x=x,y=y,color=color,size=size), data=df, show_guide=FALSE, ... )
  }
  # size
  else if( !is.null(mapping$size)){
    df$size <- get.vertex.attribute(graph,mapping$size)
    ret <- geom_point( aes(x=x,y=y,size=size), data=df, show_guide=FALSE, ... )
  }
  # color
  else if( !is.null(mapping$colour) ){
    df$colour <- get.vertex.attribute(graph,mapping$colour)
    ret <- geom_point( aes(x=x,y=y,colour=colour), data=df, show_guide=FALSE, ... )
  }
  # none 
  else
    ret <- geom_point( aes(x=x,y=y), data=df, show_guide=FALSE, ... )
  
  return( ret )

}



