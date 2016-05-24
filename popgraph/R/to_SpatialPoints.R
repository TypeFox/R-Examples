#' Translate object into a SpatialPoints object
#' 
#' Returns spatial points object
#' @param x The object containin coordinates.
#' @param stratum The name of the variable in \code{x} that represents the
#'  stratum to be used as points.
#' @param longitude The key for the attribute representing Longitude (default="Longitude")
#' @param latitude The key for the attribute representing Latitude (default="Latitude")
#' @param ... Optional arguments passed to overriden objects
#' @return An object of type SpatialPoints
#' @author Rodney J. Dyer \email{rjdyer@@vcu.edu}
#' @export
to_SpatialPoints<- function( x, stratum="Name", longitude="Longitude", latitude="Latitude", ...) {
  if( !inherits(x,"popgraph"))
    stop("This function requires a popgraph object to function")

  
  vertex.attr <- list.vertex.attributes( x )
  if( !(latitude %in% vertex.attr ) | !(longitude %in% vertex.attr) )
    stop("Your graph should have Latitude and Longitude in it before we can make it a Spatial* object.")
  
  coords <- cbind( x=get.vertex.attribute( x, longitude ),
                   y=get.vertex.attribute( x, latitude) )
  rownames( coords ) <- V(x)$name 
  pts <- SpatialPoints(coords)
  
  return( pts )
}