#' @name rain_max_world
#' @title Record rainfall totals from around the world
#' 
#' @description  World record point precipitation measurements from NOAA's 
#' National Weather Service, Hydrometeorological Design Studies Center. NOAA 
#' indicates that some records have not been verified.
#' 
#' @details
#' Variables:
#' 
#' \itemize{
#'  \item Duration - duration of the rainfall event
#'  \item Amount_in	- total rainfall amount in inches
#'  \item Amount_mm	- total rainfall amount in millimeters
#'  \item Location - region and country names
#'  \item Lat	- latitude in decimal degrees
#'  \item Lon	- longitude in decimal degrees
#'  \item Start_Date - starting date of the rainfall event in yyyy-mm-dd format	
#'  \item Estimate - Yes if the record has not been verified, blank otherwise
#' }
#'
#' @references World record point precipitation measurements, \url{http://www.nws.noaa.gov/oh/hdsc/record_precip/record_precip_world.html}, extracted May 18 2014.
#' 
#' @author Gopi Goteti
#' 
#' @docType data
#' 
#' @usage data(rain_max_world)
#' 
#' @format Data frame with 8 columns and 47 rows
#' 
#' @keywords datasets
NULL

#' @name rain_max_usa
#' @title Record rainfall totals from the USA
#' 
#' @description  USA record point precipitation measurements from NOAA's 
#' National Weather Service, Hydrometeorological Design Studies Center. NOAA 
#' indicates that some records have not been verified.
#' 
#' @details
#' Variables:
#' 
#' \itemize{
#'  \item Duration - duration of the rainfall event
#'  \item Amount_in  - total rainfall amount in inches
#'  \item Amount_mm	- total rainfall amount in millimeters
#'  \item Location - region and country names
#'  \item Lat	- latitude in decimal degrees
#'  \item Lon	- longitude in decimal degrees
#'  \item Start_Date - starting date of the rainfall event in yyyy-mm-dd format	
#'  \item Estimate - Yes or blank; Yes if the record has not been verified.
#' }
#'
#' @references USA record point precipitation measurements, \url{http://www.nws.noaa.gov/oh/hdsc/record_precip/record_precip_us.html}, extracted May 18 2014.
#' 
#' @author Gopi Goteti
#' 
#' @docType data
#' 
#' @usage data(rain_max_usa)
#' 
#' @format Data frame with 8 columns and 19 rows
#' 
#' @keywords datasets
NULL


#' @name rainfreq
#' 
#' @title Rainfall frequency estimates for the USA from the NOAA National Weather 
#' Service (NWS) division Hydrometeorological Design Studies Center (HDSC).
#'
#' @details Data from NOAA NWS is available in various formats. \pkg{rainfreq}
#' provides functionality to access the GIS format files provided by NWS' PF 
#' Data Server. \pkg{rainfreq} regional selection criterion is currently limited 
#' to the 50 states (plus DC). \pkg{rainfreq} also comes with datasets on record 
#' point rainfall measurements provided by NWS-HDSC
#' \url{http://www.nws.noaa.gov/oh/hdsc/record_precip/record_precip.html}
#'
#' @import RCurl SDMTools
#' 
#' @docType package
#' 
#' @author Gopi Goteti
#' 
#' @references NOAA NWS HDSC \url{http://www.nws.noaa.gov/oh/hdsc/index.html}. 
#' Data in GIS format is available from the PF Data Server 
#' \url{http://hdsc.nws.noaa.gov/hdsc/pfds/pfds_gis.html}.
NULL
