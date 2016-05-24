#' Prepares weather station data
#' @param WSdata             csv file with weather station data
#' @param ...                additional parameter to pass to read.csv()
#' @param height             weather station sensors height in meters
#' @param lat                latitude of weather station in decimal degrees. 
#' Negative values for south latitude
#' @param long               longitude of weather station in decimal degrees. 
#' Negative values for west longitude
#' @param elev               elevation of weather station in meters
#' @param columns            columns order of needed data. Vector containing 
#' "date", "time", "radiation", "wind", "RH" and "temp". Other values are 
#' ignored. If you have a column with date and time in the same column, you can
#' include "datetime" and "date" and "time" are no longer needed.
#' @param date.format        date format. See strptime format argument.
#' @param time.format        time format. See strptime format argument.
#' @param datetime.format    datetime format. See strptime format argument.
#' @param tz                 timezone of the weather station dates. If not 
#' present assumes the same timezone as the computer running the code. See 
#' strptime for details. 
#' @param cf                 conversion factor to convert radiation, wind, and 
#' temperature to W/m2; m/s and Celsius. See Details.
#' @param MTL                Metadata file. If not provided will look for one on
#' working directory. If provided or present will calculate weather conditions
#' on satellite flyby.
#' @details 
#' For cf, if your data is in W/m2, km/h and Celsius (radiation, wind, 
#' temperature), cf should be: cf = c(1,0.2777778,1)
#' @examples 
#' csvfile <- system.file("extdata", "apples.csv", package="water")
#' MTLfile <- system.file("extdata", "L7.MTL.txt", package="water")
#' WS <- read.WSdata(WSdata = csvfile, date.format = "%d/%m/%Y", 
#'                   lat=-35.42222, long= -71.38639, elev=201, height= 2.2,
#'                   MTL = MTLfile)
#' print(WS)
#' plot(WS, alldata=FALSE)
#' plot(WS, alldata=TRUE)
#' @author Guillermo Federico Olmedo
#' @return waterWeatherStation object, with data.frames with all data, hourly
#' data and conditions at satellite flyby.
#' @references 
#' Landsat 7 Metadata example file available from the U.S. Geological Survey.
#' Weather Station example file courtesy of CITRA, Universidad de Talca, Chile
#' @export
read.WSdata <- function(WSdata, ..., height = 2.2, lat, long, elev,
                        columns = c("date", "time", "radiation", "wind", NA,
                                    "RH", "temp", "pp"),
                        date.format = "%Y-%m-%d", time.format = "%H:%M:%S", 
                        datetime.format = "%Y-%m-%d %H:%M:%S", tz = "",
                        cf = c(1, 1, 1), MTL){
  result <- list()
  result$location <- data.frame(lat=lat, long=long, elev=elev, height=height)
  WSdata <- utils::read.csv(WSdata, ...)
  if("date" %in% columns & "time" %in% columns){
    datetime  <- paste(WSdata[, which(columns == "date")], 
                       WSdata[, which(columns == "time")])
    datetime <- strptime(datetime, format = paste(date.format, time.format), 
                         tz = tz)
  } else {
    if("datetime" %in% columns){
      datetime  <- WSdata[, which(columns == "datetime")]
      datetime <- strptime(datetime, format = datetime.format, tz = tz)
    } else {message("ERROR: date and time or datetime are needed columns")}
  }
  radiation = WSdata[, which(columns == "radiation")] * cf[1]
  wind =  WSdata[, which(columns == "wind")] * cf[2]
  RH =  WSdata[, which(columns == "RH")]
  temp =  WSdata[, which(columns == "temp")] * cf[3]
  ea = (RH/100)*0.6108*exp((17.27*temp)/(temp+237.3))
  WSdata <- data.frame(datetime=datetime, radiation=radiation, wind=wind,
                       RH=RH, ea=ea, temp=temp)
  result$alldata <- WSdata
  result$hourly <- WSdata[datetime$min==0,] 
  ## Join with satellite data
  if(missing(MTL)){MTL <- list.files(pattern = "MTL.txt", full.names = T)}
  if(length(MTL)!=0){
    MTL <- readLines(MTL, warn=FALSE)
    time.line <- grep("SCENE_CENTER_TIME",MTL,value=TRUE)
    date.line <- grep("DATE_ACQUIRED",MTL,value=TRUE)
    sat.time <-regmatches(time.line,regexec(text=time.line,
           pattern="([0-9]{2})(:)([0-9]{2})(:)([0-9]{2})(.)([0-9]{2})"))[[1]][1]
    sat.date <-regmatches(date.line,regexec(text=date.line,
                        pattern="([0-9]{4})(-)([0-9]{2})(-)([0-9]{2})"))[[1]][1]
    sat.datetime <- strptime(paste(sat.date, sat.time), 
                             format = "%Y-%m-%d %H:%M:%S", tz="GMT")
    WS.prev<-WSdata[WSdata$datetime == tail(datetime[datetime < 
                                                       sat.datetime],1),]
    WS.after <- WSdata[WSdata$datetime == datetime[datetime > 
                                                     sat.datetime][1],]
    delta1 <- as.numeric(difftime(WS.after$datetime, 
                                  WS.prev$datetime, units="secs"))
    delta2 <- as.numeric(difftime(sat.datetime, WS.prev$datetime, units="secs"))
    sat.data <- WS.prev + (WS.after - WS.prev)/delta1 * delta2
    sat.data[,2:6] <- round(sat.data[,2:6],2)
    result$at.sat <- sat.data
  }
  class(result) <- "waterWeatherStation"
  return(result)
}


#' Prepares weather station data 2
#' @param WSdata             csv file with weather station data
#' @param ...                additional parameter to pass to read.csv()
#' @param height             weather station sensors height in meters
#' @param lat                latitude of weather station in decimal degrees. 
#' Negative values for south latitude
#' @param long               longitude of weather station in decimal degrees. 
#' Negative values for west longitude
#' @param elev               elevation of weather station in meters
#' @param columns            columns order of needed data. Vector containing 
#' "date", "time", "radiation", "wind", "RH" and "temp". Other values are 
#' ignored. If you have a column with date and time in the same column, you can
#' include "datetime" and "date" and "time" are no longer needed.
#' @param date.format        date format. See strptime format argument.
#' @param time.format        time format. See strptime format argument.
#' @param datetime.format    datetime format. See strptime format argument.
#' @param tz                 timezone of the weather station dates. If not 
#' present assumes the same timezone as the computer running the code. See 
#' strptime for details. 
#' @param cf                 conversion factor to convert radiation, wind, and 
#' temperature to W/m2; m/s and Celsius. See Details.
#' @param MTL                Metadata file. If not provided will look for one on
#' working directory. If provided or present will calculate weather conditions
#' on satellite flyby.
#' @details 
#' For cf, if your data is in W/m2, km/h and Celsius (radiation, wind, 
#' temperature), cf should be: cf = c(1,0.2777778,1)
#' @author Guillermo Federico Olmedo
#' @export
read.WSdata2 <- function(WSdata, ..., height = 2.2, lat, long, elev,
                         columns = c("date", "time", "radiation", "wind", NA,
                                     "RH", "temp", "pp"),
                         date.format = "%d/%m/%Y", time.format = "%H:%M:%S", 
                         datetime.format = "%Y-%m-%d %H:%M:%S", tz = "",
                         cf = c(1, 3.6, 1, 1), MTL){
  read.WSdata(WSdata, ..., height = 2.2, lat = lat, long = long, elev = elev,
              columns = columns, date.format = date.format, 
              time.format = time.format, datetime.format = datetime.format,
              tz = tz, cf = cf, MTL)
}


#' Plot method for waterWeatherStation S3 class
#' @param x       waterWeatherStation object. See read.WSdata()
#' @param hourly  If TRUE will plot only hourly data, instead of all data 
#' from the waterWeatherStation object
#' @param sat     If TRUE, and if the waterWeatherStation object was created
#' using a Landsat Metadata File, will plot only data from day of the satellite
#' overpass 
#' @param ...      additional parameters to pass to plot()
#' @author Guillermo Federico Olmedo
#' @importFrom utils read.csv
#' @importFrom graphics abline axis axis.POSIXct mtext par points
#' @export
#' @method plot waterWeatherStation
plot.waterWeatherStation <- function(x, hourly=FALSE, sat=TRUE, ...){
  # Based on http://evolvingspaces.blogspot.cl/2011/05/multiple-y-axis-in-r-plot.html
  WSp <- x$hourly
  atsat  <- as.POSIXct(x$at.sat$datetime)
  if(hourly == FALSE) {WSp <- x$alldata}
  if(sat == TRUE & atsat > 0){WSp <- WSp[as.Date(WSp$datetime) == as.Date(atsat),]}
  time <- WSp$datetime
  graphics::par(mar=c(5, 7, 4, 7) + 0.1)
  plot(time, WSp$radiation, axes=F, ylim=c(0,max(WSp$radiation)), xlab="", 
       ylab="",type="l",col="red", main="",xlim=range(time), ...)
  graphics::abline(v=atsat, lwd=5, col="gray")
  graphics::text(atsat, max(WSp$radiation)*0.85, "Satellite Flyby", cex=0.7, 
                 adj=c(NA, -0.5), srt=90, col="gray")
  graphics::points(time,WSp$radiation,pch=20,col="red")
  graphics::axis(2, ylim=c(0,max(WSp$radiation)),col="red",lwd=1, cex.axis=0.5)
  graphics::mtext(2,text="Solar radiation (W.m-2)",line=1.7, cex=0.7)
  graphics::par(new=T)
  plot(time, WSp$wind, axes=F, ylim=c(0,max(WSp$wind)), xlab="", ylab="", 
       type="l",col="green",lty=2, main="",xlim=range(time),lwd=2,...)
  graphics::axis(2, ylim=c(0,max(WSp$wind)),col="green",lwd=1,line=3.5, cex.axis=0.5)
  graphics::points(time, WSp$wind,pch=20,col="green")
  graphics::mtext(2,text="Wind speed (m.s-1)",line=5.5, cex=0.7)
  graphics::par(new=T)
  plot(time, WSp$temp, axes=F, ylim=c(0,max(WSp$temp)), xlab="", ylab="", 
       type="l",col="black",lty=2, main="",xlim=range(time),lwd=2,...)
  graphics::axis(4, ylim=c(0,max(WSp$temp)),col="black",lwd=1, cex.axis=0.5)
  graphics::points(time, WSp$temp,pch=20,col="black")
  graphics::mtext(4,text="Temperature (C)",line=1.7, cex=0.7)
  graphics::par(new=T)
  plot(time, WSp$ea, axes=F, ylim=c(0,max(WSp$ea)), xlab="", ylab="", 
       type="l",col="blue",lty=2, main="",xlim=range(time),lwd=2, ...)
  graphics::axis(4, ylim=c(0,max(WSp$ea)),col="blue",lwd=1,line=3.5, cex.axis=0.5)
  graphics::points(time, WSp$ea,pch=20,col="blue")
  graphics::mtext(4,text="vapor pressure (kPa)",line=5.5, cex=0.7)
  graphics::axis.POSIXct(1, range(time))
  graphics::mtext("Time",side=1,col="black",line=2)
}

#' Print method for waterWeatherStation S3 class
#' @param x        waterWeatherStation object. See read.WSdata()
#' @param ...      additional parameters to pass to print()    
#' @author Guillermo Federico Olmedo
#' @export
#' @method print waterWeatherStation
print.waterWeatherStation <- function(x, ...){
  cat("Weather Station at lat:", round(x$location$lat, 2), "long:", 
      round(x$location$long, 2), "elev:", round(x$location$elev, 2), "\n")
  cat("Summary:\n")
  print(summary(x$alldata[,2:6]), ...)
  if(!is.null(x$at.sat)){
  cat("\n Conditions at satellite flyby:\n")
  print(x$at.sat)}
}


#' Export data.frame from waterWeatherStation Object
#' @description 
#' Export weather conditions at satellite flyby from waterWeatherStation object
#' @param WeatherStation    waterWeatherStation object. See read.WSdata()
#' @author Guillermo Federico Olmedo
#' @export
getDataWS <- function(WeatherStation){
  date <- as.POSIXlt(WeatherStation$at.sat$datetime, format="%Y-%m-%d %H:%M:%S")
  WeatherStation2 <- data.frame(wind=WeatherStation$at.sat$wind, 
                               RH=WeatherStation$at.sat$RH, 
                               temp=WeatherStation$at.sat$temp, 
                               radiation=WeatherStation$at.sat$radiation, 
                               height=WeatherStation$location$height,
                               lat=WeatherStation$location$lat, 
                               long=WeatherStation$location$long, 
                               elev=WeatherStation$location$elev, 
                               DOY=date$yday+1, hours=date$hour + date$min/60 + 
                                 date$sec/3600)
  return(WeatherStation2)
}


