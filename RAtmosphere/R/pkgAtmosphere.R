standard_atmosphere<-function(z,vout=3,T0=288.15,P0=101325,R0=1.225){
  #function developed with information provided from:
  # http://scipp.ucsc.edu/outreach/balloon/atmos/1976 Standard Atmosphere.htm
  #
  # z is altitude over sea level
  #
  # parameter vout has default value 3
  #  1 for temperature, 2 for pressure, 3 for density

  level0<-c(0,11,20,32,47,51) *1000.
  level1<-c(11,20,32,47,51,71)*1000.

  func<-vector("list",length=3)#ncol=3,nrow=3)
  func[[1]]<-vector("list",length=6)#ncol=3,nrow=3)
  func[[2]]<-vector("list",length=6)#ncol=3,nrow=3)
  func[[3]]<-vector("list",length=6)

  func[[1]][[1]]<-function(z,T0=288.15){return(T0*(1.-z/44329))}
  func[[1]][[2]]<-function(z,T0=288.15){return(T0*(0.751865))}
  func[[1]][[3]]<-function(z,T0=288.15){return(T0*(0.682457 + z/288136))}
  func[[1]][[4]]<-function(z,T0=288.15){return(T0*(0.482561 + z/102906))}
  func[[1]][[5]]<-function(z,T0=288.15){return(T0*(0.939268))}
  func[[1]][[6]]<-function(z,T0=288.15){return(T0*(1.434843 - z/102906))}
  func[[2]][[1]]<-function(z,P0=101325){return(P0*(1. - z/44329)^5.255876)}
  func[[2]][[2]]<-function(z,P0=101325){return(P0*(0.223361)*exp((10999-z)/6341.4))}
  func[[2]][[3]]<-function(z,P0=101325){return(P0*(0.988626 + z/198903)^(-34.16319))}
  func[[2]][[4]]<-function(z,P0=101325){return(P0*(0.898309 + z/55280)^(-12.20114))}
  func[[2]][[5]]<-function(z,P0=101325){return(P0*(0.00109456)*exp((46998-z)/7922))}
  func[[2]][[6]]<-function(z,P0=101325){return(P0*(0.838263 - z/176142)^(12.20114))}
  func[[3]][[1]]<-function(z,R0=1.225){return(R0*(1. - z/44329)^4.255876)}
  func[[3]][[2]]<-function(z,R0=1.225){return(R0*(0.297076)*exp((10999-z)/6341.4))}
  func[[3]][[3]]<-function(z,R0=1.225){return(R0*(0.978261 + z/201010)^(-35.16319))}
  func[[3]][[4]]<-function(z,R0=1.225){return(R0*(0.857003 + z/57944)^(-13.20114))}
  func[[3]][[5]]<-function(z,R0=1.225){return(R0*(0.00116533)*exp((46998-z)/7922 ))}
  func[[3]][[6]]<-function(z,R0=1.225){return(R0*(0.79899 - z/184800)^11.20114)}

  out<-z
  for(i in 1:length(z)){
   n_level<-min(which(level1>=z[i]),na.rm=TRUE)
#   print(c(vout,n_level))
   out[i]<-func[[vout]][[n_level]](z[i])
  }
  return(out)
}


betamol<-function(P,T,lambda=1064.)2.938*10^4.1053*P/T/lambda^4.0117

betamol_standard<-function(z,lambda=1064.,T0=288.15,P0=101325.,R0=1.225){
betamol(standard_atmosphere(z,vout=2,T0=T0,P0=P0,R0=R0)/100.,standard_atmosphere(z,vout=1,T0=T0,P0=P0,R0=R0),lambda=lambda)

}

SZA<-function(timein=Sys.time(),Lat = 50.910335,Lon = 11.568740){
# Calculate solar zenith angle
# according to http://solardat.uoregon.edu/SolarRadiationBasics.html
# calculations have been modified for positive East longitudes and time in UTC
# Extract time information
#( hour, minute, second, dummy, n ) = time.gmtime()[3:8]
# Calculate declination of the sun d
#if (is.vector(timein)){

 sza<-vector("numeric",length=length(timein))
  for(i in 1:length(timein)){
    time<-as.POSIXlt(timein[i],tz='GMT')
d2r=pi/180.
r2d=1./d2r
d<-23.45 * d2r * sin(d2r *360. *(284. + time$yday) / 365.) # [rad]
#print (sprintf("d = %f", d * 180 / pi))
# Calculate equation of time
if (time$yday <= 106){
    E_qt <- -14.2 * sin(pi * (time$yday + 7.) / 111.)      # Eq. SR.4a [minutes]
}else{if (time$yday<= 166){
    E_qt <-   4.0 * sin(pi * (time$yday - 106.) / 59.)     # Eq. SR.4b [minutes]
}else{if (time$yday<= 246){
    E_qt <-  -6.5 * sin(pi * (time$yday - 166.) / 80.)     # Eq. SR.4c [minutes]
}else{E_qt <-  16.4 * sin(pi * (time$yday - 247.) / 113.)}}}    # Eq. SR.4d [minutes]
# Get UTC time T
T<-time$hour + time$min / 60.0 + time$sec / 3600.0 # [hours]
# Calculate solar time T_solar (East longitudes are positive!)
Longitude<-Lon#*d2r
T_solar<-T + Longitude / 15. + E_qt / 60. # [hours]
# Calculate hour angle w (positive: from midnight to noon, negative: from noon to midnight)
w<-pi * (12. - T_solar) / 12. # [rad]
# Calculate solar zenith angle Z
l<-Lat * d2r # [rad]

sza[i]<-90.-asin(sin(l) * sin(d) + cos(l) * cos(d) * cos(w)) * r2d # [deg]

  }
#}
return(sza)
}

suncalc<-function(d,Lat=0.,Long=0.,UTC=TRUE){
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
  rad<-function(x)pi*x/180.

  ##Radius of the earth (km)
  R=6378.

  ##Radians between the xy-plane and the ecliptic plane
  epsilon=rad(23.45)

  ##Convert observer's latitude to radians
  L=rad(Lat)

  ## Calculate offset of sunrise based on longitude (min)
  ## If Long is negative, then the mod represents degrees West of
  ## a standard time meridian, so timing of sunrise and sunset should
  ## be made later.
  if (UTC){
           timezone = 0
	   }else{
		   timezone =-4*(abs(Long)%%15)*sign(Long)
	   }

  ## The earth's mean distance from the sun (km)
  r = 149598000.

  theta = 2*pi/365.25*(d-80.)

  z.s = r*sin(theta)*sin(epsilon)
  r.p = sqrt(r^2-z.s^2)

  t0 = 1440./(2.*pi)*acos((R-z.s*sin(L))/(r.p*cos(L)))

  ##a kludge adjustment for the radius of the sun
  that = t0+5.

  ## Adjust "noon" for the fact that the earth's orbit is not circular:
  n = 720.-10.*sin(4*pi*(d-80)/365.25)+8.*sin(2.*pi*d/365.25)

  ## now sunrise and sunset are:
  if(UTC){
  sunrise = (n-that)/60.
  sunset = (n+that)/60.
}else{
  sunrise = (n-that+timezone)/60.
  sunset = (n+that+timezone)/60.
}

  return(list("sunrise" = sunrise,"sunset" = sunset))
}
