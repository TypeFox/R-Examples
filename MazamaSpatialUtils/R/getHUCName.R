#' @keywords locator
#' @export
#' @title Return HUC Names at Specified Locations
#' @param lon vector of longitudes in decimal degrees
#' @param lat vector of latitudes in decimal degrees
#' @param dataset name of spatial dataset to use
#' @param HUCs vector of Hydrologic Unit Codes
#' @param allData logical specifying whether to return a full dataframe
#' @description Uses spatial comparison to determine which HUC polygons the 
#'     locations fall into and returns the HUC names for those polygons.
#'     
#'     If \code{allData=TRUE}, additional data is returned.
#' @return Vector of HUC names. 
#' @seealso getSpatialData
 

getHUCName <- function(lon, lat, dataset='WBDHU10-ms', HUCs=NULL, allData=FALSE) {
  
  # Sanity check
  if (!exists(dataset)) {
    stop('Missing database. Please loadSpatialData("',dataset,'")',call.=FALSE)
  }
  
  # Use standard internal name (assumes pre-loaded dataset)
  SPDF <- get(dataset) 
  
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
    
    return(HUCName)
    
  }
  
  
}

