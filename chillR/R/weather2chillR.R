#' Convert downloaded weather to chillR format
#' 
#' Convert downloaded weather data into a data frame that makes running other
#' chillR functions easy.
#' 
#' weather databases, from which chillR can download data: NOAA NCDC Global
#' Summary of the Day - "GSOD"
#' (https://data.noaa.gov/dataset/global-surface-summary-of-the-day-gsod)
#' 
#' Weather Underground ("Wunderground") airport database - "Wunderground"
#' (http://www.wunderground.com/)
#' 
#' California Irrigation Management Information System (CIMIS) - "CIMIS"
#' (http://www.cimis.water.ca.gov/)
#' 
#' University of California Integrated Pest Management (UCIPM) - "UCIPM"
#' (http://ipm.ucdavis.edu/WEATHER/)
#' 
#' data should first be downloaded with get_weather. Then the database name is
#' passed to the function and can be skipped in the call. If only a data.frame
#' is provided, then the database name must be specified.
#' 
#' Processing the data with this function will make the data work well with the
#' remainder of this package.
#' 
#' @param downloaded_weather weather file downloaded with the get_weather
#' function. This can be a data.frame or a list with elements database and
#' weather as produced by get_weather
#' @param database weather database that the file was downloaded from. Can only
#' be "GSOD" at this point.
#' @param drop_most boolean variable indicating if most columns should be
#' dropped from the file. If set to TRUE (default), only essential columns for
#' running chillR functions are retained.
#' @return a data.frame with weather data, according to the downloaded file
#' provided as input. If drop_most is FALSE, all columns from the original
#' dataset are preserved, although some column names are adjusted to chillR's
#' preferences ("Year","Month","Day","Tmin","Tmax","Tmean","Prec"). If
#' drop_most is TRUE, only columns likely to be of interest to chillR users are
#' retained. If a list with elements database and weather is passed to this
#' function, this structure will be retained in the output.
#' @note Many databases have data quality flags, which may sometimes indicate
#' that data aren't reliable. These are not considered by this function!
#' @author Eike Luedeling
#' @references The chillR package:
#' 
#' Luedeling E, Kunz A and Blanke M, 2013. Identification of chilling and heat
#' requirements of cherry trees - a statistical approach. International Journal
#' of Biometeorology 57,679-689.
#' @keywords utilities
#' @examples
#' 
#' #stat_list<-get_weather(location=c(lat=40,lon=-120,ele=150),time_interval=c(2015,2016),
#' #database="UCIPM")
#' #chillRcode<-stat_list[which(stat_list$Perc_interval_covered==
#' #max(stat_list$Perc_interval_covered)),"chillR_code"][1]
#' #chillRcode should equal "DOYLE.C" now.
#' gw<-get_weather(location="DOYLE.C",time_interval=c(2002,2002),database="UCIPM")
#' #weather<-weather2chillR(gw$weather,"GSOD")
#' weather<-weather2chillR(gw)
#' 
#' @export weather2chillR
weather2chillR<-function(downloaded_weather,database="GSOD",drop_most=TRUE)
{dw<-downloaded_weather
if(is.list(dw)) if(names(dw)[1]=="database") database=dw$database

if(database=="GSOD")  
  return(handle_gsod(dw,drop_most=drop_most))
 
  if(database=="CIMIS")
      return(handle_cimis(dw,drop_most=drop_most))

if(database=="Wunderground")
  return(handle_wunderground(dw,drop_most=drop_most))

if(database=="UCIPM")
  return(handle_ucipm(dw,drop_most=drop_most))


}
