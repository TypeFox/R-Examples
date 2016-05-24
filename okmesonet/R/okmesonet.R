## roxygen documentation for okmesonet package and subsequent datasets

#' Oklahoma Mesonet climatological data
#'
#' The \href{http://www.mesonet.org/}{Oklahoma Mesonet} is a network of 
#' automated climate monitoring stations throughout Oklahoma, USA. 
#' Data collection began January 01, 1994; as of November 2014, there are 120 
#' active stations. Measurements are recorded every five minutes and sent to a 
#' central facility for verification and quality control by the Oklahoma 
#' Climatological Survey.
#'
#' The list of Oklahoma Mesonet stations was updated November 28, 2014. 
#' Use \code{\link{updatestn}} to manually update the list.
#'
#' Data access may be restricted by organization and/or location. Please refer 
#' to and follow policies found within the
#' \href{http://www.mesonet.org/index.php/site/about/data_access_and_pricing}{Oklahoma 
#' Mesonet Data Access Policy}. The authors and maintainers of \code{okmesonet}
#' assume no responsibility for the use or misuse of \code{okmesonet}.
#'
#' @name okmesonet-package
#' @aliases okmesonet
#' @docType package
NULL
#' Bison locations
#'
#' Location information for one American bison (\emph{Bison bison}) at The
#' Nature Conservancy Tallgrass Prairie Preserve, Jan 31 - Feb 5, 2011.
#' Information includes time, latitude, and longitude.
#' @docType data
#' @name bison
#' @format Data frame
NULL
#' List of Oklahoma Mesonet stations
#'
#' List of \href{http://www.mesonet.org/}{Oklahoma Mesonet} stations. 
#' The Oklahoma Mesonet is a network of 
#' automated climate monitoring stations throughout Oklahoma, USA. 
#' Data collection began January 01, 1994; as of November 2014, there are 120 
#' active stations. Measurements are recorded every five minutes and sent to a 
#' central facility for verification and quality control by the Oklahoma 
#' Climatological Survey.
#'
#' Station list is updated each time \code{okmesonet} is loaded. 
#'
#' Variables associated with each station include:
#' \itemize{
#'  \item Identifier: four letter station identifier
#'  \item Number: station identification number
#'  \item Name: station name
#'  \item Town: closest town
#'  \item County: Oklahoma county
#'  \item Latitude: station latitude, decimal degrees
#'  \item Longitude: station longitude, decimal degrees
#'  \item Elevation: station elevation, meters
#'  \item Commissioned: date station was commissioned
#'  \item Decommissioned: date station was decommissioned. "2099-12-31" 
#'    represents currently active
#' }
#'
#' @docType data
#' @name okstations
#' @format Data frame
#' @seealso \code{\link{updatestn}} to manually update station list.
NULL
#' Check geomeso.csv
#'
#' Internal function to check access to geomeso.csv
#' @docType data
#' @name checkgeoinfo
NULL
#' Download Oklahoma Mesonet station list
#'
#' Internal function to download Oklahoma Mesonet station list
#' @docType data
#' @name downloadstn
NULL

checkgeoinfo <- function() {
  ## Check to verify internet connectivity and access to 
  ## http://www.mesonet.org/index.php/api/siteinfo/from_all_active_with_geo_fields/format/csv/ 
  ## is available.
  ##
  ## Arguments: none
  ## Returns: logical indicating success
  path.geoinfo <- paste("http://www.mesonet.org/index.php/api/siteinfo/",
                        "from_all_active_with_geo_fields/format/csv/", sep ="")
  internetoption <- getOption("internet.info")
  options(show.error.messages=FALSE, internet.info=3)
  checkconn <- try(readLines(path.geoinfo, n=1))
  options(show.error.messages=TRUE, internet.info=internetoption)
  if (inherits(checkconn,"try-error")) return(FALSE) else return(TRUE)
}

downloadstn <- function() {
  ## Retrieve Oklahoma Mesonet station information, including:
  ## station identifier, station number, station name, nearest town, latitude,
  ## longitude, elevation, date commissioned, date decommissioned
  ##
  ## Arguments: none
  ## Returns: data frame with above information
  
  ## use .checkgeoinfo() to check connectivity
  path.geoinfo <- paste("http://www.mesonet.org/index.php/api/siteinfo/",
                        "from_all_active_with_geo_fields/format/csv/", sep ="")
  if(checkgeoinfo()==FALSE) {
    warn.msg <- paste("Oklahoma Mesonet station list unavailable. Check", 
                      path.geoinfo, "for", 
                      "connectivity and run updatestn() to update station", 
                      "list.")
    warning(warn.msg, call.=F)
  } else {
    geoinfo <- read.csv(path.geoinfo, header=T, as.is=T)
    stationlist <- data.frame(Identifier=geoinfo$stid, Number=geoinfo$stnm,
                              Name=geoinfo$name, Town=geoinfo$city, 
                              County=geoinfo$cnty, Latitude=geoinfo$nlat, 
                              Longitude=geoinfo$elon, Elevation=geoinfo$elev,
                              Commissioned=strptime(geoinfo$datc, 
                                                    format="%Y%m%d",
                                                    tz="America/Chicago"),
                              Decommissioned=strptime(geoinfo$datd, 
                                                      format="%Y%m%d",
                                                      tz="America/Chicago"),
                              stringsAsFactors=F)
    return(stationlist)
  }
}