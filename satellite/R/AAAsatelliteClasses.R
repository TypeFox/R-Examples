#' An S4 class to represent a satellite data file
#'
#' @slot name name of the data file without extension
#' @slot filepath full path and file of the data file
#' @slot path path to the data file
#' @slot file filename incl. extension of the data file
#' @slot extension extension of the data file
#' 
#' @exportClass SatelliteInfo
#' 
setClass("SatelliteInfo", 
         representation(
           name ="character",
           filepath = "character",
           path = "character",
           file = "character",
           extension = "character"
         )
)


#' An S4 class to represent satellite metadata
#'
#' @slot meta a data frame object containing the data
#' 
#' @exportClass SatelliteMetaData
#' 
setClass("SatelliteMetaData",
         representation(
           meta = "data.frame"
         )
)


#' An S4 class to represent satellite data
#'
#' @slot layers a list object containing individual RasterLayer objects
#' 
#' @exportClass SatelliteLayers
#' 
setClass("SatelliteLayers",
         representation(
           layers = "list"
         )
)


#' An S4 class to represent satellite log data
#'
#' @slot log a list object containing information on individual processing steps
#' 
#' @exportClass SatelliteLog
#' 
setClass("SatelliteLog", 
         representation(
           log = "list"
         )
)

#' An S4 class to represent a complete satellite dataset
#' 
#' @exportClass Satellite
#' 
setClass("Satellite", 
         contains = c("SatelliteLayers", "SatelliteMetaData", "SatelliteLog")
)         
