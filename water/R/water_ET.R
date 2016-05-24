#' Calculates ET using Penman Monteith hourly formula
#' @param WeatherStation     a data frame with all the needed fields (see 
#' example)
#' @param hours              time of the day in hours in 24hs format
#' @param DOY                day of year
#' @param long.z             longitude for local time
#' @param ET.instantaneous   Logical. True if you want to calculate 
#' instantaneous ET instead of hourly ET. See Details.
#' @param ET                 "ETo" or "ETr"
#' @param height             weather station sensors height in meters
#' @param lat                latitude of weather station in decimal degrees. 
#' Negative values for south latitude
#' @param long               longitude of weather station in decimal degrees. 
#' Negative values for west longitude
#' @param elev               elevation of weather station in meters
#' @details 
#' The only difference on instantaneous ET is how the hour is interpreted. On 
#' FALSE, and for example at 11:00, ET is calculated between 10:00 and 11:00, 
#' on TRUE Et is calculated at 11:00 hs.
#' @return ET hourly in mm.h-1
#' @author Guillermo Federico Olmedo
#' @examples 
#' WeatherStation  <- data.frame(wind=4.72,
#'                               RH=59, 
#'                               temp=24.3,
#'                               radiation=675, 
#'                               height=2.2, 
#'                               lat=-35.37, 
#'                               long=71.5946, 
#'                               elev=124)
#'   hourlyET(WeatherStation, hours=10.5, DOY=363, long.z=71.635)
#' @export
#' @references 
#' Allen 2005 ASCE
hourlyET <- function(WeatherStation, hours, DOY, long.z=WeatherStation$long, 
                     ET.instantaneous=FALSE, ET="ETo", height=2, lat, long, elev){
  if(class(WeatherStation)== "waterWeatherStation"){
    if(!is.null(WeatherStation$at.sat)){
      WeatherStation <- getDataWS(WeatherStation)
      ET.instantaneous = TRUE
      DOY <- WeatherStation$DOY
      hours <-  WeatherStation$hours}
  } else {
    if("datetime" %in% names(WeatherStation)){
      date <- as.POSIXlt(WeatherStation$datetime, format="%Y-%m-%d %H:%M:%S")
      DOY=date$yday+1
      hours=date$hour + date$min/60 + date$sec/3600}
    if(!missing(lat)){WeatherStation$lat <- lat}
    if(!missing(long)){WeatherStation$long <- long}
    if(!missing(elev)){WeatherStation$elev <- elev}
    if(!missing(height)){WeatherStation$height <- height}
    if(length(WeatherStation$DOY) != 0 & length(WeatherStation$hours) != 0 &
       missing(DOY) & missing(hours)){
      DOY <- WeatherStation$DOY
      hours <- WeatherStation$hours
    } 
  }
  tempK <- WeatherStation$temp + 273.16
  Rs <- WeatherStation$radiation * 3600 / 1e6
  P <- 101.3*((293-0.0065*WeatherStation$elev)/293)^5.26
  psi <- 0.000665*P
  Delta <- 2503 * exp((17.27*WeatherStation$temp)/
                  (WeatherStation$temp+237.3))/((WeatherStation$temp+237.3)^2)
  ea.sat <- 0.61078*exp((17.269*WeatherStation$temp)/
                          (WeatherStation$temp+237.3))
  ea <- (WeatherStation$RH/100)*ea.sat
  DPV <- ea.sat - ea
  dr <- 1 + 0.033*(cos(2*pi*DOY/365))
  delta <- 0.409*sin((2*pi*DOY/365)-1.39)
  phi <- WeatherStation$lat*(pi/180)
  b <- 2*pi*(DOY-81)/364
  Sc <- 0.1645*sin(2*b)-0.1255*cos(b)-0.025*sin(b)
  if(ET.instantaneous==TRUE){hours <- hours +0.5}
  hour.angle <- (pi/12)*((hours-0.5+0.06667*(WeatherStation$long*
                                               pi/180-long.z*pi/180)+Sc)-12)
  w1 <- hour.angle-((pi)/24)
  w2 <- hour.angle+((pi)/24)
  hour.angle.s <- acos(-tan(phi)*tan(delta))
  w1c <- ifelse(w1< -hour.angle.s, -hour.angle.s, 
                ifelse(w1>hour.angle.s, hour.angle.s, ifelse(w1>w2, w2, w1)))
  w2c <- ifelse(w2< -hour.angle.s, -hour.angle.s, 
                ifelse(w2>hour.angle.s, hour.angle.s, w2))
  Beta <- asin((sin(phi)*sin(delta)+cos(phi)*cos(delta)*cos(hour.angle)))
  Ra <- ifelse(Beta <= 0, 1e-45, ((12/pi)*4.92*dr)*
                 (((w2c-w1c)*sin(phi)*sin(delta))+(cos(phi)*cos(delta)*
                                                     (sin(w2)-sin(w1)))))
  Rso <- (0.75+2e-5*WeatherStation$elev)*Ra
  Rs.Rso <- ifelse(Rs/Rso<=0.3, 0, ifelse(Rs/Rso>=1, 1, Rs/Rso))
  fcd <- ifelse(1.35*Rs.Rso-0.35<=0.05, 0.05, 
                ifelse(1.35*Rs.Rso-0.35<1, 1.35*Rs.Rso-0.35,1))
  Rn.a <- ((1-0.23)*Rs) - (2.042e-10*fcd*(0.34-0.14*(ea^0.5))*tempK^4)
  Ghr <- c(0.1, 0.5)
  if(ET=="ETr"){Ghr <- c(0.04, 0.2)}
  G.day <- Rn.a * ifelse(Rn.a > 0, Ghr[1], Ghr[2])
  wind.2 <- WeatherStation$wind *(4.87/(log(67.8*WeatherStation$height-5.42)))
  windfactors <- c(0.24, 0.96)
  if(ET=="ETr"){windfactors <-  c(0.25, 1.7)}
  windf <- ifelse(WeatherStation$radiation > 0, windfactors[1], windfactors[2])
  ET.hourly <- ((0.408*Delta*(Rn.a-G.day))+(psi*(37/tempK)*wind.2*(DPV)))/
    (Delta+(psi*(1+(windf*wind.2))))
  if(ET=="ETr"){
    ET.hourly <- ((0.408*Delta*(Rn.a-G.day))+(psi*(66/tempK)*wind.2*(DPV)))/
      (Delta+(psi*(1+(windf*wind.2))))
  }
  return(ET.hourly)
}

#' Calculates ET-24hs from energy balance and Weather Station 
#' @param Rn             Net radiation. See netRadiation()
#' @param G              Soil Heat Flux. See soilHeatFlux()
#' @param H              Sensible Heat Flux. See calcH()
#' @param Ts             Land surface temperature. See surfaceTemperature()
#' @param WeatherStation WeatherStation data at the flyby from the satellite. 
#' Can be a waterWeatherStation object calculate using read.WSdata and MTL file
#' @param ETr.daily      hourly ETr for every hour of the day. See dailyET()
#' @param C.rad          correction term used in sloping terrain to correct for 
#' variation in 24 h versus instantaneous energy availability. See Allen (2007)
#' @author Guillermo Federico Olmedo
#' @references 
#' R. G. Allen, M. Tasumi, and R. Trezza, "Satellite-based energy balance for mapping evapotranspiration with internalized calibration (METRIC) - Model" Journal of Irrigation and Drainage Engineering, vol. 133, p. 380, 2007
#' @importFrom grDevices colorRampPalette
#' @export
ET24h <- function(Rn, G, H, Ts, WeatherStation, ETr.daily, C.rad=1){
  if(class(WeatherStation)== "waterWeatherStation"){
    WeatherStation <- getDataWS(WeatherStation)
  }
  LE = Rn - G - H
  ET.inst <- 3600*LE/((2.501 - 0.00236 * (Ts - 273.15)) * (1e6))
  ETo.hourly <- hourlyET(WeatherStation, WeatherStation$hours, WeatherStation$DOY)
  ETr.Fr <- ET.inst/ETo.hourly
  ET.24 <- ETr.Fr * ETr.daily * C.rad
  #ET.24[ET.24 < 0]  <- 0
  #ET.24[ET.24 > quantile(ET.24, 0.9)] <- quantile(ET.24, 0.9)
  rgb.palette <- grDevices::colorRampPalette(c("red3","snow2","blue"),  space = "rgb")
  print(spplot(ET.24, col.regions=rgb.palette, main= "24-Hour Evapotranspiration (mm/day)",
               colorkey=list(height=1), at=seq(0,ceiling(ETr.daily*1.5),length.out=50), maxpixels=ncell(ET.24) * 0.3))
  saveLoadClean(imagestack = ET.24, 
                file = "ET24", overwrite=TRUE)
}


#' Calculates daily ET using Penman Monteith hourly formula for every hour
#' @param WeatherStation a data frame with all the needed fields (see example)
#' @param DOY      day of year
#' @param height   weather station sensors height in meters
#' @param lat      latitude in decimal degrees of the weather station
#' @param long     longitude in decimal degrees of the weather station
#' @param elev     elevation in meters of the weather station
#' @param ET       "ETo" or "ETr"
#' @param long.z   longitude for local time
#' @return ET      daily in mm.h-1
#' @author Guillermo Federico Olmedo
#' @export
#' @references 
#' Allen 2005 ASCE
#' @examples 
#' csvfile <- system.file("extdata", "apples.csv", package="water")
#' 
#' WeatherStation <- read.WSdata(WSdata = csvfile, date.format = "%d/%m/%Y", 
#' lat=-35.42222, long= -71.38639, elev=201, height= 2.2, cf=c(1,0.2777778,1,1))
#' 
#' dailyET(WeatherStation = WeatherStation, lat=-35.422, long=-71.386, elev=124, 
#' ET="ETo")
#' 
dailyET <- function(WeatherStation, DOY, height, lat, long, elev, ET="ETo", 
                       long.z=WeatherStation$long){
  if(class(WeatherStation)== "waterWeatherStation"){
    if(missing(height)){height <- WeatherStation$location$height}
    if(missing(lat)){lat <- WeatherStation$location$lat}
    if(missing(long)){long <- WeatherStation$location$long}
    if(missing(elev)){elev <- WeatherStation$location$elev}
    ET.daily <- vector()
    for(i in 1:24){
      date <- as.POSIXlt(WeatherStation$hourly[i,1], format="%Y-%m-%d %H:%M:%S")
      ET.daily <- c(ET.daily, hourlyET(WeatherStation$hourly[i,], lat=lat, 
                                      long = long, elev=elev, ET=ET, 
                                      height = height))
      }
  } else {
  print("not yet")}
  return(sum(ET.daily))
}
