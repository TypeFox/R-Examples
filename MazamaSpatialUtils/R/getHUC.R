#' @keywords locator
#' @export
#' @title Return HUCs at Specified Locations
#' @param lon vector of longitudes in decimal degrees
#' @param lat vector of latitudes in decimal degrees
#' @param SPDF spatial polygons dataset of HUCs
#' @param HUCs vector of Hydrologic Unit Codes
#' @param allData logical specifying whether to return a full dataframe
#' @description Uses spatial comparison to determine which HUC polygons the 
#'     locations fall into and returns the HUC identifier strings for those polygons.
#'     
#'     If \code{allData=TRUE}, additional data is returned.
#' @return Vector of HUC identifiers.
#' @seealso getSpatialData
 

getHUC <- function(lon, lat, SPDF, HUCs=NULL, allData=FALSE) {
  
  # TODO:  sanity check that dataset is a character string 
  
#   # Sanity check
#   if (!exists(dataset)) {
#     stop('Missing database. Please loadSpatialData("',dataset,'")',call.=FALSE)
#   }
#   
#   # Use standard internal name (assumes pre-loaded dataset)
#   SPDF <- get(dataset) 
#   
  # Identify HUC string partial matches to use as a mask 
  if (!is.null(HUCs)){
    HUCMask <- rep(FALSE, nrow(SPDF))
    for (HUC in HUCs){
      regex <- paste0('^', HUC)
      mask <- stringr::str_detect(SPDF@data$HUC, regex)
      HUCMask <- HUCMask | mask
    }
    SPDF <- SPDF[HUCMask,]
  }
  
  # Pull out rows from SPDF@data based on whether a lon-lat point falls into a certain polygon 
  locationsDF <- getSpatialData(lon,lat,SPDF)
  
  if (allData) {
    
    return(locationsDF)
    
  } else {
    
    HUC <- locationsDF$HUC
    HUCName <- locationsDF$HUCName
    
    return(HUC)
    
  }
  
  
}

