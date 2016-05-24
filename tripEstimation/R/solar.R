"solar" <-
function(day) {
  ## day - times as POSIXct
  
  ## Extract components of time (GMT)
  tm <- as.POSIXlt(day,tz="GMT")           
  hh <- tm$hour
  mm <- tm$min
  ss <- tm$sec

  ## Time as Julian day
  jday <- julday(tm)+(hh+(mm+ss/60)/60)/24
  
  ## Time as Julian century
  t <- (jday-2451545)/36525

  ## Geometric mean anomaly for the sun (degrees)
  M <- 357.52911+t*(35999.05029-0.0001537*t)

  ## Equation of centre for the sun (degrees)
  eqcent <- sin(pi/180*M)*(1.914602-t*(0.004817+0.000014*t))+
    sin(pi/180*2*M)*(0.019993-0.000101*t)+
      sin(pi/180*3*M)*0.000289
    
  ## The geometric mean sun longitude (degrees)
  L0 <- 280.46646+t*(36000.76983+0.0003032*t)
  ## Limit to [0,360)
  L0 <- L0%%360 

  ## The true longitude of the sun (degrees)
  lambda0 <- L0 + eqcent
  
  ## The apparent longitude of the sun (degrees)
  omega <- 125.04-1934.136*t
  lambda <- lambda0-0.00569-0.00478*sin(pi/180*omega)
  

  ## The mean obliquity of the ecliptic (degrees)
  seconds <- 21.448-t*(46.815+t*(0.00059-t*(0.001813)))
  obliq0 <- 23+(26+(seconds/60))/60 

  ## The corrected obliquity of the ecliptic (degrees)
  omega <- 125.04-1934.136*t
  obliq <- obliq0 + 0.00256*cos(pi/180*omega)

  ## The eccentricity of earth's orbit
  e <- 0.016708634-t*(0.000042037+0.0000001267*t)
  
  ## The equation of time (minutes of time)
  y <- tan(pi/180*obliq/2)^2
  eqtime <- 180/pi*4*(y*sin(pi/180*2*L0) -
                      2*e*sin(pi/180*M) +
                      4*e*y*sin(pi/180*M)*cos(pi/180*2*L0) -
                      0.5*y^2*sin(pi/180*4*L0) -
                      1.25*e^2*sin(pi/180*2*M))

  ## The sun's declination (radians)
  solarDec <- asin(sin(pi/180*obliq)*sin(pi/180*lambda))
  sinSolarDec <- sin(solarDec)
  cosSolarDec <- cos(solarDec)

  ## ALT
  ## sinSolarDec <- sin(pi/180*obliq)*sin(pi/180*lambda)
  ## cosSolarDec <- cos(asin(sinSolarDec))

  ## Solar time unadjusted for longitude (degrees)
  solarTime <- (hh*60+mm+ss/60+eqtime)/4

  ## Return solar constants
  list(solarTime=solarTime,
       sinSolarDec=sinSolarDec,
       cosSolarDec=cosSolarDec)
}

