#' A simple way to access maps stored in the package.
#' 
#' A simple way to access maps stored in the package.
#' 
#' 
#' @param resolution options "coarse","low","less islands","li","high". For
#' "high" you need to install the package rworldxtra
#' @param projection DEPRECATED OCTOBER 2012 to reproject maps see spTransform
#' in rgdal
#' @return A SpatialPolygonsDataFrame object.
#' @author Barry Rowlingson & Andy South
#' @keywords misc
#' @examples
#' 
#' plot(getMap())
#' 
#' @export getMap
`getMap` <-

function(resolution="coarse",projection=NA){
  
  
  resolutionOptions <- c("coarse","low","less islands","li","high") 
  
  if( ! resolution %in% resolutionOptions)
  {
    warning("resolution should be set to one of :",paste(resolutionOptions,""),"\nsetting to coarse as default\n")
    resolution="coarse"
  }
    
  
  if (! is.na(projection)) 
     {
      warning("the projection argument to getMap() in rworldmap is deprecated and will be removed in a future release.
			        Returning an unprojected map, use spTransform() from package rgdal to project.")
     }

  #coarsest resolution map
 if (resolution == "coarse") {
    data("countriesCoarse", envir = environment(),package = "rworldmap")
    mapWithData <- get("countriesCoarse", envir = environment(), inherits=FALSE) 
    
  }  else if (resolution == "low") {
    data("countriesLow", envir = environment(),package = "rworldmap")
    mapWithData <- get("countriesLow", envir = environment(), inherits=FALSE) 
    
  }  else if (resolution == "less islands" | resolution == "li") {
    data("countriesCoarseLessIslands", envir = environment(),package = "rworldmap")
    mapWithData <- get("countriesCoarseLessIslands", envir = environment(), inherits=FALSE)   
    
  }  else if (resolution == "high" && !requireNamespace(package = "rworldxtra", quietly=TRUE)) {
    warning("for resolution='high' option you need to install package rworldxtra, using low resolution version for now")
    data("countriesLow", envir = environment(), package = "rworldmap")
    mapWithData <- get("countriesLow", envir = environment(), inherits=FALSE)
    
  }  else if (resolution == "high" ) {
    data("countriesHigh", envir = environment(), package = "rworldxtra")
    mapWithData <- get("countriesHigh", envir = environment(), inherits=FALSE)    
  }  

 
  #trying eez map : does work - temporarily removed 31/8 to get permission
  #if (resolution == "eez") {
  #  data("eezMap", envir = environment(), package = "rworldmap")
  #  mapWithData <- get("eezMap")
  #} 
  #else
  
  return(mapWithData)
}
