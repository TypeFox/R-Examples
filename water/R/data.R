#' Landsat 7 scene from east of Talca, Chile
#'
#' A RasterStack object with a Landsat 7 ETM+ subset image. 
#' Band names are stored in layer names. Images are USGS L1 raw data. 
#'
#' Metadata for this image is also provided as system file. You can
#' load it with:  system.file("extdata", "L7.MTL.txt", package="water")
#' 
#' Data available from the U.S. Geological Survey.
#'
#' @source \url{http://www.usgs.gov/}
"L7_Talca"


#' SRTM DEM from east of Talca, Chile
#'
#' A RasterLayer object with a a Digital Elevation Model. 
#' 
#' Original data is float. Example data is integer.
#' 
#' Data available from the U.S. Geological Survey.
#'
#' @source \url{http://www.usgs.gov/}
"DEM_Talca"