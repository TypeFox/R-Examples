## Functions copied from SGAT and temporariliy included into GeoLight:
## 1) solar
## 2) zenith
## 3) refracted
## 4) twilight
## 5) geolight.convert



##' Calculate solar time, the equation of time and solar declination
##'
##' The solar time, the equation of time and the sine and cosine of
##' the solar declination are calculted for the times specified by
##' \code{tm} using the same methods as
##' \url{www.esrl.noaa.gov/gmd/grad/solcalc/}.
##' @title Solar Time and Declination
##' @param tm a vector of POSIXct times.
##' @return A list containing the following vectors.
##' \item{\code{solarTime}}{the solar time (degrees)}
##' \item{\code{eqnTime}}{the equation of time (minutes of time)}
##' \item{\code{sinSolarDec}}{sine of the solar declination}
##' \item{\code{cosSolarDec}}{cosine of the solar declination}
##' @seealso \code{\link{zenith}}
##' @examples
##' ## Current solar time
##' solar(Sys.time())
##' @export
solar <- function(tm) {
  
  rad <- pi/180
  
  ## Time as Julian day (R form)
  Jd <- as.numeric(tm)/86400.0+2440587.5
  
  ## Time as Julian century [G]
  Jc <- (Jd-2451545)/36525
  
  ## The geometric mean sun longitude (degrees) [I]
  L0 <- (280.46646+Jc*(36000.76983+0.0003032*Jc))%%360
  
  
  ## Geometric mean anomaly for the sun (degrees) [J]
  M <- 357.52911+Jc*(35999.05029-0.0001537*Jc)
  
  ## The eccentricity of earth's orbit [K]
  e <- 0.016708634-Jc*(0.000042037+0.0000001267*Jc)
  
  ## Equation of centre for the sun (degrees) [L]
  eqctr <- sin(rad*M)*(1.914602-Jc*(0.004817+0.000014*Jc))+
    sin(rad*2*M)*(0.019993-0.000101*Jc)+
    sin(rad*3*M)*0.000289
  
  ## The true longitude of the sun (degrees) [M]
  lambda0 <- L0 + eqctr
  
  ## The apparent longitude of the sun (degrees) [P]
  omega <- 125.04-1934.136*Jc
  lambda <- lambda0-0.00569-0.00478*sin(rad*omega)
  
  
  ## The mean obliquity of the ecliptic (degrees) [Q]
  seconds <- 21.448-Jc*(46.815+Jc*(0.00059-Jc*(0.001813)))
  obliq0 <- 23+(26+(seconds/60))/60
  
  ## The corrected obliquity of the ecliptic (degrees) [R]
  omega <- 125.04-1934.136*Jc
  obliq <- obliq0 + 0.00256*cos(rad*omega)
  
  ## The equation of time (minutes of time) [U,V]
  y <- tan(rad*obliq/2)^2
  eqnTime <- 4/rad*(y*sin(rad*2*L0) -
                      2*e*sin(rad*M) +
                      4*e*y*sin(rad*M)*cos(rad*2*L0) -
                      0.5*y^2*sin(rad*4*L0) -
                      1.25*e^2*sin(rad*2*M))
  
  ## The sun's declination (radians) [T]
  solarDec <- asin(sin(rad*obliq)*sin(rad*lambda))
  sinSolarDec <- sin(solarDec)
  cosSolarDec <- cos(solarDec)
  
  ## Solar time unadjusted for longitude (degrees) [AB!!]
  ## Am missing a mod 360 here, but is only used within cosine.
  solarTime <- ((Jd-0.5)%%1*1440+eqnTime)/4
  #solarTime <- ((Jd-2440587.5)*1440+eqnTime)/4
  
  ## Return solar constants
  list(solarTime=solarTime,
       eqnTime=eqnTime,
       sinSolarDec=sinSolarDec,
       cosSolarDec=cosSolarDec)
}



##' Calculate the solar zenith angle for given times and locations
##'
##' \code{zenith} uses the solar time and declination calculated by
##' \code{solar} to compute the solar zenith angle for given times and
##' locations, using the same methods as
##' \url{www.esrl.noaa.gov/gmd/grad/solcalc/}.  This function does not
##' adjust for atmospheric refraction see \code{\link{refracted}}.
##' @title Solar Zenith Angle
##' @param sun list of solar time and declination computed by \code{solar}.
##' @param lon vector of longitudes.
##' @param lat vector latitudes.
##' @return A vector of solar zenith angles (degrees) for the given
##' locations and times.
##' @seealso \code{\link{solar}}
##' @examples
##' ## Approx location of Sydney Harbour Bridge
##' lon <- 151.211
##' lat <- -33.852
##' ## Solar zenith angle for noon on the first of May 2000
##' ## at the Sydney Harbour Bridge
##' s <- solar(as.POSIXct("2000-05-01 12:00:00","EST"))
##' zenith(s,lon,lat)
##' @export
zenith <- function(sun,lon,lat) {
  
  rad <- pi/180
  
  ## Suns hour angle (degrees) [AC!!]
  hourAngle <- sun$solarTime+lon-180
  #hourAngle <- sun$solarTime%%360+lon-180
  
  ## Cosine of sun's zenith [AD]
  cosZenith <- (sin(rad*lat)*sun$sinSolarDec+
                  cos(rad*lat)*sun$cosSolarDec*cos(rad*hourAngle))
  
  ## Limit to [-1,1] [!!]
  cosZenith[cosZenith > 1] <- 1
  cosZenith[cosZenith < -1] <- -1
  
  ## Ignore refraction correction
  acos(cosZenith)/rad
}


##' Adjust the solar zenith angle for atmospheric refraction.
##'
##' Given a vector of solar zeniths computed by \code{\link{zenith}},
##' \code{refracted} calculates the solar zeniths adjusted for the
##' effect of atmospheric refraction.
##'
##' \code{unrefracted} is the inverse of \code{refracted}. Given a
##' (single) solar zenith adjusted for the effect of atmospheric
##' refraction, \code{unrefracted} calculates the solar zenith as
##' computed by \code{\link{zenith}}.
##'
##' @title Atmospheric Refraction
##' @param zenith zenith angle (degrees) to adjust.
##' @return vector of zenith angles (degrees) adjusted for atmospheric
##' refraction.
##' @examples
##' ## Refraction causes the sun to appears higher on the horizon
##' refracted(85:92)
##' ## unrefracted gives unadjusted zenith (see SGAT)
##' 
##' @export
refracted <- function(zenith) {
  rad <- pi/180
  elev <- 90-zenith
  te <- tan((rad)*elev)
  ## Atmospheric Refraction [AF]
  r <- ifelse(elev>85,0,
              ifelse(elev>5,58.1/te-0.07/te^3+0.000086/te^5,
                     ifelse(elev>-0.575,
                            1735+elev*(-518.2+elev*(103.4+elev*(-12.79+elev*0.711))),-20.772/te)))
  ## Corrected Zenith [90-AG]
  zenith-r/3600
}





##' Estimate time of sunrsie or sunset for a given day and location
##'
##' \code{twilight} uses an iterative algorithm to estimate times of
##' sunrise and sunset.
##'
##' Note that these functions return the twilight that occurs on the
##' same date GMT as \code{tm}, and so sunset may occur before
##' sunrise, depending upon latitude.
##'
##' Solar declination and equation of time vary slowly over the day,
##' and so the values of the Solar declination and equation of time at
##' sunrise/sunset are well approximated by their values at 6AM/6PM
##' local time. The sun's hour angle and hence sunrise/sunset for the
##' required zenith can then be caclulates from these approximations.
##' The calculation is then repeated using the approximate
##' sunrise/sunset times to derive more accurate values of the Solar
##' declination and equation of time and hence better approximations
##' of sunrise/sunset.  The process is repreated and is accurate to
##' less than 2 seconds within 2 or 3 iterations.
##'
##' \code{sunrise} and \code{sunset} are simple wrappers for \code{twilight}.
##' @title Times of Sunrise and Sunset
##' @param tm vector of approximate times of twilight.
##' @param lon vector of longitudes.
##' @param lat vector of latitudes.
##' @param rise logical vector indicating whether to compute rise or set.
##' @param zenith the solar zenith angle that defines twilight.
##' @param iters number of iteratve refinements made to the initial
##' approximation.
##' @return a vector of twilight times.
##' @export
twilight <- function(tm,lon,lat,rise,zenith=96,iters=3) {
  
  ## Compute date
  date <- as.POSIXlt(tm)
  date$hour <- date$min <- date$sec <- 0
  date <- as.POSIXct(date,"GMT")
  
  lon <- (lon+180)%%360-180
  ## GMT equivalent of 6am or 6pm local time
  twl <- date+240*(ifelse(rise,90,270)-lon)
  ## Iteratively improve estimate
  for(k in seq_len(iters)) {
    s <- solar(twl)
    s$solarTime <- s$solarTime%%360
    solarTime <- 4*twilight.solartime(s,lon,lat,rise,zenith)-s$eqnTime
    twl <- date+60*solarTime
  }
  twl
}


twilight.solartime <- function(solar,lon,lat,rise,zenith=96) {
  rad <- pi/180
  cosz <- cos(rad*zenith)
  cosHA <- (cosz-sin(rad*lat)*solar$sinSolarDec)/(cos(rad*lat)*solar$cosSolarDec)
  ## Compute the sun's hour angle from its declination for this location
  hourAngle <- ifelse(rise,360,0)+ifelse(rise,-1,1)*suppressWarnings(acos(cosHA)/rad)
  ## Solar time of sunrise at this zenith angle, lon and lat
  #(hourAngle+180-lon)%%360
  #360*(solar$solarTime%/%360)+solarTime
  solarTime <- (hourAngle+180-lon)%%360
  (solarTime-solar$solarTime+180)%%360-180+solar$solarTime
}


i.geolight.convert <- function(tFirst,tSecond,type) {
  tm <- .POSIXct(c(as.POSIXct(tFirst,"GMT"),
                   as.POSIXct(tSecond,"GMT")),"GMT")
  keep <- !duplicated(tm)
  tm <- tm[keep]
  rise <- c(type==1,type!=1)[keep]
  ord <- order(tm)
  data.frame(Twilight=tm[ord],Rise=rise[ord])
}


##' Convert GeoLight data
##'
##' This function converts from the tFirst, tSecond format used by
##' GeoLight to the twilight, rise format used by Stella and Estelle.
##' @title Convert GeoLight Format
##' @param tFirst times of first twilight.
##' @param tSecond times of second twilight.
##' @param type type of twilight.
##' @return A data frame with columns
##' \item{\code{twilight}}{times of twilight as POSIXct objects.}
##' \item{\code{rise}}{logical vector indicating which twilights are sunrise.}
##' @export
geolight.convert <- function(tFirst,tSecond,type) {
  tm <- .POSIXct(c(as.POSIXct(tFirst,"GMT"),
                   as.POSIXct(tSecond,"GMT")),"GMT")
  keep <- !duplicated(tm)
  tm <- tm[keep]
  rise <- c(type==1,type!=1)[keep]
  ord <- order(tm)
  data.frame(Twilight=tm[ord],Rise=rise[ord])
}

