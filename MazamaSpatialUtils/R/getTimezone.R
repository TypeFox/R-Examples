#' @keywords locator
#' @export
#' @title Return Olson Timezones at Specified Locations
#' @param lon vector of longitudes in decimal degrees
#' @param lat vector of latitudes in decimal degrees
#' @param dataset name of spatial dataset to use
#' @param countryCodes vector of countryCodes
#' @param allData logical specifying whether to return a full dataframe
#' @param useBuffering logical flag specyfing the use of location buffering to find the nearest polygon if not target polygon is found
#' @description Uses spatial comparison to determine which timezone polygons the 
#'     locations fall into and returns the Olson timezone strings for those polygons.
#'     
#'     Specification of \code{countryCodes} limits spatial searching to the specified
#'     countries and greatly improves performance.
#'     
#'     If \code{allData=TRUE}, additional data is returned.
#' @return Vector of Olson timezones.
#' @examples
#' lon <- seq(-120,-60,5)
#' lat <- seq(20,80,5)
#' getTimezone(lon,lat)
#' @references \url{http://efele.net/maps/tz/}
#' @seealso SimpleTimezones
#' @seealso getSpatialData
getTimezone <- function(lon, lat, dataset="SimpleTimezones", countryCodes=NULL, allData=FALSE, useBuffering=FALSE) {
  
  # Sanity check
  if (!exists(dataset)) {
    stop('Missing database. Please loadSpatialData("',dataset,'")',call.=FALSE)
  }
 
  SPDF <- get(dataset)
  
  # Subset by country before searching
  if (!is.null(countryCodes)) SPDF <- SPDF[SPDF$countryCode %in% countryCodes,]
  
  SPDF <- getSpatialData(lon,lat,SPDF,useBuffering=useBuffering)
  
  if (allData) {

    return(SPDF)
    
  } else {
    
    timezone <- SPDF$timezone
    
    # Sanity check -- missing timezone implies location over water  
    badMask <- is.na(timezone)
    if (sum(badMask) > 0) {
      if(is.null(countryCodes)) {
        warning(paste(sum(badMask),"locations appear to be over international waters and no timezone can be assigned"))
      } else {
        warning(paste(sum(badMask),"locations appear to be either over international waters or not in given countryCodes and no timezone can be assigned"))
      }
    }  
    
    return(timezone)
    
  }
  
}

