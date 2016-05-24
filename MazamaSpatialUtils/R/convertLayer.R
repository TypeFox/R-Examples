#' @keywords internal
#' @export
#' @title Convert Shapefile Layer to Spatial Polygon Dataframe
#' @param dsn dsn argument to readOGR
#' @param layerName layer argument to readOGR
#' @param encoding encoding string (.e.g. 'latin1') passed to rgdal::readOGR()
#' @description Raw shapefiles are read in using the \code{rgdal::readOGR()} function from the \pkg{rgdal} package.
#' Spatial data are reprojected onto a standard projection with \code{"+proj=longlat +ellps=GRS80 +datum=NAD83 +no_defs"} before being returned.
#' @return An object of class \code{SpatialPolygonsDataFrame}
convertLayer <- function(dsn="", layerName="", encoding=NULL) {
  
  # readOGR does not interpret '~' so do that with dirname()
  dsn <- path.expand(dsn)
  
  # Use package internal data directory
  dataDir <- getSpatialDataDir()
  
  # Always require a data directory
  if (is.null(dataDir)) {
    stop('dataDir must be specified and must be a user writable directory', call.=FALSE)
  }
  
  # Switch directories
  oldDir <- getwd()
  setwd(dataDir)
  
  # Load the shapefiles 
  data_projected <- rgdal::readOGR(dsn=dsn, layer=layerName, stringsAsFactors=FALSE, encoding=encoding)
  
  # Return to user directory
  setwd(oldDir)
  
  # Reproject to standard projection
  SPDF <- sp::spTransform(data_projected, sp::CRS("+proj=longlat +ellps=GRS80 +datum=NAD83 +no_defs"))
  
  return(SPDF)  
}

