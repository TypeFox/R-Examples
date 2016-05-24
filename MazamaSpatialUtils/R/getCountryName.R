#' @keywords locator
#' @export
#' @title Return Country Names at Specified Locations
#' @param lon vector of longitudes in decimal degrees
#' @param lat vector of latitudes in decimal degrees
#' @param dataset name of spatial dataset to use
#' @param countryCodes vector of countryCodes
#' @param allData logical specifying whether a full dataframe should be returned
#' @param useBuffering logical flag specyfing the use of location buffering to find the nearest polygon if not target polygon is found
#' @description Uses spatial comparison to determine which country polygons the 
#'     locations fall into and returns the country name for those polygons.
#'     
#'     If \code{allData=TRUE}, additional data is returned.
#' @return Vector of country names in English.
#' @examples
#' lon <- seq(0,50)
#' lat <- seq(0,50)
#' getCountryName(lon,lat)
#' @references \url{http://www.naturalearthdata.com/downloads/10m-cultural-vectors/}
#' @seealso SimpleCountries
#' @seealso getSpatialData
getCountryName <- function(lon, lat, dataset='SimpleCountries', countryCodes=NULL, allData=FALSE, useBuffering=FALSE) {
  
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
    
    countryName <- locationsDF$countryName
    
    # Sanity check -- missing countryName implies location over water  
    badMask <- is.na(countryName)
    if (sum(badMask) > 0) {
      if(is.null(countryCodes)) {
        warning(paste(sum(badMask),"locations appear to be over international waters and no country can be assigned"))
      } else {
        warning(paste(sum(badMask),"locations appear to be either over international waters or not in given countryCodes and no country can be assigned"))
      }
    }  
    
    return(countryName)
    
  }
  
  
}

