#' Function to create a lakeMorpho class - this is input to all other methods
#' 
#' This is a helper function that creates a lakeMorpho class object
#' 
#' @param inLake input lake SpatialPolygons object
#' @param inElev input elevation model raster object
#' @param inCatch input catchement SpatialPolygons object, can be buffer
#'        around lake
#' @param inLakeDist input euclidean distance raster that measures distance
#'        from shore to any pixel in the lake
#' @param lakeOnEdge Boolean indicating if inCatch (or lake Buffer) extends 
#'        beyond extent of elevation data
#'          
#' @export
#' @return Returns an objec of class 'lakeMorpho' 
#' @seealso lakeSurroundTopo    


# May need to be done as a method (i.e. no need to @export) TO DO: Add null place holders for all possible
# lakeMorpho metrics (eg various lines)
lakeMorphoClass <- function(inLake, inElev, inCatch, inLakeDist, lakeOnEdge = F) {
    lmorpho <- list(lake = inLake, elev = inElev, surround = inCatch, lakeDistance = inLakeDist, lakeOnEdge = lakeOnEdge)
    class(lmorpho) <- "lakeMorpho"
    return(lmorpho)
} 
