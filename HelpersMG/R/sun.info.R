#' sun.info estimate the time of sunrise and sunset according to longitude, latitude and date
#' @title Estimate the time of sunrise and sunset according to longitude, latitude and date
#' @author Marc Girondot \email{marc.girondot@@u-psud.fr}
#' @return A data.frame with information about daily sun
#' @param latitude The latitude at which estimate the sun fates
#' @param longitude The longitude at which estimate the sun fates
#' @param date  A vector with the time at which sun fates are needed
#' @family Periodic patterns of indices
#' @description Estimate the sun fates according to latitude and date.\cr
#' Based on Teets, D.A. 2003. Predicting sunrise and sunset times. The College Mathematics Journal 34(4):317-321.\cr
#' Can be compared with the function \code{sunrise.set()} of package \code{StreamMetabolism}.
#' @examples
#' \dontrun{
#' # Generate a timeserie of time
#' date <- seq(from=as.Date("2000-01-01"), to=as.Date("2000-12-31"), by="1 day")
#' plot(date, sun.info(date, latitude=23, longitude=0)$day.length, bty="n", 
#'  las=1, type="l", xlab="Ordinal days", ylab="Day length in hours")
#' plot(date, sun.info(date, latitude=23, longitude=0)$sunrise, bty="n", 
#'  las=1, type="l", xlab="Ordinal days", ylab="Sun rise in hours")
#' }
#' @export

sun.info <- function(date, latitude, longitude){
  
  d <- as.numeric(as.POSIXlt(date)$yday)+1
  Lat <- latitude
  Long <- longitude
  
  ## d is the day of year
  ## Lat is latitude in decimal degrees
  ## Long is longitude in decimal degrees (negative == West)
  
  ##This method is copied from:
  ##Teets, D.A. 2003. Predicting sunrise and sunset times.
  ##  The College Mathematics Journal 34(4):317-321.
  
  ## At the default location the estimates of sunrise and sunset are within
  ## seven minutes of the correct times (http://aa.usno.navy.mil/data/docs/RS_OneYear.php)
  ## with a mean of 2.4 minutes error.
  
  ## Function to convert degrees to radians
  rad <- function(x) pi*x/180
  
  ##Radius of the earth (km)
  R=6378
  
  ##Radians between the xy-plane and the ecliptic plane
  epsilon=rad(23.45)
  
  ##Convert observer's latitude to radians
  L=rad(Lat)
  
  ## Calculate offset of sunrise based on longitude (min)
  ## If Long is negative, then the mod represents degrees West of
  ## a standard time meridian, so timing of sunrise and sunset should
  ## be made later.
  timezone = -4*(abs(Long)%%15)*sign(Long)
  
  ## The earth's mean distance from the sun (km)
  r = 149598000
  
  theta = 2*pi/365.25*(d-80)
  
  z.s = r*sin(theta)*sin(epsilon)
  r.p = sqrt(r^2-z.s^2)
  
  t0 = 1440/(2*pi)*acos((R-z.s*sin(L))/(r.p*cos(L)))
  
  ##a kludge adjustment for the radius of the sun
  that = t0+5 
  
  ## Adjust "noon" for the fact that the earth's orbit is not circular:
  n = 720-10*sin(4*pi*(d-80)/365.25)+8*sin(2*pi*d/365.25)
  
  ## now sunrise and sunset are:
  sunrise = (n-that+timezone)/60
  sunset = (n+that+timezone)/60
  
  UTC <- (((7.5+Long)%%360)) %/% 15
  if (UTC>12) {UTC <- 12-UTC;tz <- "Etc/GMT"} else {tz <- "Etc/GMT+"}
  tz <- paste0(tz, UTC)
  
  df <- data.frame(sunrise = sunrise, sunset = sunset, day.length=sunset- sunrise, 
                   date.time.sunrise=as.POSIXlt(format(date, "%Y-%m-%d"), tz=tz)+sunrise*60*60, 
                   date.time.sunset=as.POSIXlt(format(date, "%Y-%m-%d"), tz=tz)+sunset*60*60)
  
  sunrise.UTC <- as.POSIXlt(format(df$date.time.sunrise, format="%Y-%m-%d %H:%M:%S"), 
                           tz="UTC", use.tz=TRUE)
  sunrise.UTC.dec <- sunrise.UTC$hour + sunrise.UTC$min/60 + sunrise.UTC$sec/3600
  sunset.UTC <- as.POSIXlt(format(df$date.time.sunset, format="%Y-%m-%d %H:%M:%S"), 
                            tz="UTC", use.tz=TRUE)
  sunset.UTC.dec <- sunset.UTC$hour + sunset.UTC$min/60 + sunset.UTC$sec/3600
  df <- cbind(df, date.time.sunrise.UTC=sunrise.UTC, date.time.sunset.UTC=sunset.UTC, 
              time.sunrise.UTC=sunrise.UTC.dec, time.sunset.UTC=sunset.UTC.dec)
  
  return(df)
}

