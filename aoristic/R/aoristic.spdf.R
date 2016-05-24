#' creating a spatial.polygon.data.frame for aoristic analysis
#' @param data data.frame with a minimum of 4 columns representing FromDateTime, ToDateTime, lon, lat 
#' @param DateTimeFrom a character vector of the column name for FromDateTime (POSIXct date-time object)
#' @param DateTimeTo a character vector of the column name for ToDateTime (POSIXct date-time object).  If ending date-time is missing, the duration of an event will be coded as 1 hour.
#' @param lon a character vector longitude
#' @param lat a character vector of the column name for latitude
#' @return spatial point data frame (SPDF)
#' @import lubridate classInt reshape2 GISTools ggplot2 spatstat
#' @importFrom sp SpatialPointsDataFrame
#' @export
#' @references Ratcliffe, J. H. (2002). Aoristic Signatures and the Spatio-Temporal Analysis of High Volume Crime Patterns. Journal of Quantitative Criminology, 18(1), 23-43. 
#' @examples
#' \donttest{
#' data(aoristic)
#' data.spdf <- aoristic.spdf(data=arlington, 
#'    DateTimeFrom="DateTimeFrom", DateTimeTo="DateTimeTo", 
#'    lon="lon", lat="lat")
#' }
aoristic.spdf <- function(data, DateTimeFrom, DateTimeTo, lon, lat){
  
  # check arguments
  if(!is.data.frame(data)) {stop("the input data frame specified is not a data.frame object")}
  if(!class(data[,DateTimeFrom])[1]=="POSIXct") {stop("the DateTimeFrom field is not POSIXct object.  Use the lubridate package before using this function")}
  if(!class(data[,DateTimeTo])[1]=="POSIXct")   {stop("the DateTimeTo field is not POSIXct object.  Use the lubridate package before using this function")}
  if(!is.numeric(data[,lon])){stop("the longitude is not numeric")}
  if(!is.numeric(data[,lat])){stop("the latitude is not numeric")}
  
  # ver 0.3
  CRS <- "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"
  # ver 0.5
  # CRS <- "+init=epsg:4326" # causes errors in aoristic.density with GE_SpatialGrid()
  
  
  duration <- as.numeric(difftime(data[,DateTimeTo], data[,DateTimeFrom], units="hours") + 1 )
  HourFrom <- hour(data[,DateTimeFrom])
  
  duration[duration>24] <- 24
  duration <- ceiling(duration)
  duration[is.na(duration)] <- 1 # recode duration as 1 hour if timeTo is missing
  
  for (i in 0:23){
    assign(x=paste("time", i, sep=""), 0)
  }
  
  # create df for aoristic --------------------
  id <- seq(1,nrow(data), 1)
  
  temp.df <- data.frame(time0=numeric(),  time1=numeric(),  time2=numeric(),  time3=numeric(), time4=numeric(), time5=numeric(), time6=numeric(), time7=numeric(),
                        time8=numeric(),  time9=numeric(),  time10=numeric(), time11=numeric(), time12=numeric(), time13=numeric(), time14=numeric(), time15=numeric(),
                        time16=numeric(), time17=numeric(), time18=numeric(), time19=numeric(), time20=numeric(), time21=numeric(), time22=numeric(), time23=numeric())
  
  for (i in 1:nrow(data)){
    h <- HourFrom[i]
    d <- duration[i]
    
    if(d >= 24){ 
      temp.df[i,1:24] <- 1/d
    } else if (h+d < 24){ 
      temp.df[i,(h+1):(h+d)] <- 1/d
    } else {
      temp.df[i,(h+1):24] <- 1/d
      temp.df[i,1:(h+d-23)] <- 1/d
    }
  }
  
  temp.df[is.na(temp.df)] <- 0
  
  lon <- data[, lon]
  lat <- data[, lat]
  
  data2 <- as.data.frame(cbind(id, HourFrom, duration, temp.df, lon, lat))
  
  # create SPDF data----------
  data.spdf <- SpatialPointsDataFrame(data=data2, coords=matrix(c(data2$lon, data2$lat), ncol=2), 
                                      proj4string=CRS(CRS))
  return(data.spdf)
  
}