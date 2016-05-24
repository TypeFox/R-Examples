#' Converts edge set to SpatialLines object
#' 
#' This is a convienence function that takes the edge set and returns a 
#'  SpatialLines object.
#' @param graph An object of type \code{popgraph}. This graph must already
#'  be decroated with latitude and longitude attributes.  
#' @param latitude The name of the Latitude attribute (default="Latitude")
#' @param longitude The name of the Longitude attribute (default="Longitude")
#' @param ... Ignored
#' @return A \code{SpatialLines} object
#' @author Rodney J. Dyer <rjdyer@@vcu.edu>
#' @export
to_SpatialLines <- function( graph, latitude="Latitude", longitude="Longitude", ...) {
  if( !inherits(graph,"popgraph"))
    stop("This function requires a popgraph object to function")
  
  names <- V(graph)$name
  
  vertex.attr <- list.vertex.attributes( graph )
  if( !(latitude %in% vertex.attr ) | !(longitude %in% vertex.attr) )
    stop("Your graph should have Latitude and Longitude in it before we can make it a Spatial* object.")
  
  lat <- get.vertex.attribute( graph, latitude )
  lon <- get.vertex.attribute( graph, longitude )
  
  edgeList <- list()
  all.edges <- get.edgelist( graph )
  
  for( i in 1:nrow(all.edges) ){
    name1 <- all.edges[i,1]
    name2 <- all.edges[i,2]
    idx1 <- which( names==name1 )
    idx2 <- which( names==name2 )
  
    coord <- matrix( c( lon[idx1], lat[idx1], lon[idx2], lat[idx2] ), nrow=2, byrow=TRUE )
    colnames(coord) <- c("x","y")
    rownames( coord ) <- c( names[idx1], names[idx2] )
    edgeName <- paste( "Edge", names[idx1], names[idx2] )    
    edgeList[edgeName] <- Lines( list( Line( coord )), ID=edgeName )
  }
  return( SpatialLines( edgeList ) )
}