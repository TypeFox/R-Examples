#' @keywords locator
#' @export
#' @title Return State ISO Codes at Specified Locations
#' @param lon vector of longitudes in decimal degrees
#' @param lat vector of latitudes in decimal degrees
#' @param dataset name of spatial dataset to use
#' @param countryCodes vector of country codes
#' @param allData logical specifying whether to return a full dataframe
#' @param useBuffering logical flag specyfing the use of location buffering to find the nearest polygon if not target polygon is found
#' @description Uses spatial comparison to determine which 'state' polygons the 
#'     locations fall into and returns the ISO 3166 2-character state code
#'     strings for those polygons.
#'     
#'     Specification of \code{countryCodes} limits spatial searching to the specified
#'     countries and greatly improves performance.
#'     
#'     If \code{allData=TRUE}, additional data is returned.
#' @return Vector of ISO-3166-2 alpha-2 state codes.
#' @examples
#' \dontrun{
#' lon <- seq(-140,-90)
#' lat <- seq(20,70)
#' getStateCode(lon,lat)
#' }
#' @seealso getSpatialData
getStateCode <- function(lon, lat, dataset='NaturalEarthAdm1', countryCodes=NULL, allData=FALSE, useBuffering=FALSE) {
  
  # Sanity check
  if (!exists(dataset)) {
    stop('Missing database. Please loadSpatialData("',dataset,'")',call.=FALSE)
  }
  
  SPDF <- get(dataset)
  
  # Subset by country before searching
  if (!is.null(countryCodes)) SPDF <- SPDF[SPDF$countryCode %in% countryCodes,]
  
  locationsDF <- getSpatialData(lon,lat,SPDF,useBuffering=useBuffering)
  
  if (allData) {
    
    return(locationsDF)
    
  } else {
    
    stateCode <- locationsDF$stateCode
    
    # Sanity check -- missing stateCode implies location over water  
    badMask <- is.na(stateCode)
    if (sum(badMask) > 0) {
      if(is.null(countryCodes)) {
        warning(paste(sum(badMask),"locations appear to be over international waters and no stateCode can be assigned"))
      } else {
        warning(paste(sum(badMask),"locations appear to be either over international waters or not in given countryCodes and no stateCode can be assigned"))
      }
    }  
    
    return(stateCode)
    
  }
  
  
}

