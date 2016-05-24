#'Calculates position of the sun
#'
#'@description Calculates the zenith and azimuth angle of the position of the sun, based
#'mostly on routines from Iqbal (1983).
#'

#'
#'By default, it is assumed that the time of day is not given in local apparent
#'time (LAT, also known as 'solar time'). To convert the standard time to LAT,
#'the longitude of the location, and the longitude of the nearest time zone
#'border must be given.
#'
#'Alternatively, use \code{LAT=TRUE} to specify that the time of day is in LAT
#'(that is, solar maximum occurs exactly at noon).
#'
#'The user can specify a number of timesteps (\code{KHRS}), so that the solar
#'positions are calculated for the midpoint of each timestep (this is used
#'within YplantQMC). Alternatively, specify \code{timeofday} directly.
#'
#'@param year YYYY - to account for eccentricity (small effect).
#'@param month Month number
#'@param day Day of month number
#'@param lat Latitude, degrees.
#'@param long Longitude, degrees. Optional (only used for local apparent time
#'correction.)
#'@param tzlong Longitude of the nearest timezone border
#'@param KHRS Number of timesteps in a day (optional). See Details.
#'@param timeofday Optional, time of day (in hours) (a vector of any length) to
#'calculate the position of the sun
#'@param LAT Logical (default=FALSE). Are the times of day given in 'local
#'apparent time'?
#'@return A list with the following components: \describe{
#'\item{list("hour")}{Time in decimal hours} \item{list("altitude")}{Solar
#'altitude (degrees)} \item{list("azimuth")}{Solar azimuth (degrees. N=0,
#'E=90)} \item{list("daylength")}{Day length in hours}
#'\item{list("sunset")}{Time of sunset (hours)} \item{list("zenrad")}{Solar
#'zenith position (radians)} }
#'@note This routine is no doubt less accurate that the NOAA routines provided
#'by the \code{solarpos} function in the \code{maptools} package. It is easier
#'to use, though.
#'@author Remko Duursma, based mostly on original FORTRAN code by Belinda
#'Medlyn.
#'@seealso \code{\link{setHemi}}
#'@references Iqbal, B., 1983. An Introduction to Solar Radiation. Academic
#'Press, New York, 386 pp
#'@keywords misc
#'@examples
#'
#'
#'# Simple use
#'zenaz(month=8, day=16, timeofday=12, lat=-33)
#'
#'# Get half-hourly solar positions
#'hourpos <- zenaz(month=2, day=16, KHRS=48, lat=-33, long=155, tzlong=150)
#' with(hourpos, plot(hour, altitude, type='o',ylab=expression(Altitude~(degree))))
#'
#'@export
zenaz <- function(year=2012, month=4, day=1, 
  lat= -33.6, long=150.7, 
	tzlong=long, KHRS=24, timeofday=NA, LAT=FALSE){    

  # Private function SUN (only gets declination, time corrections.)
      SUN <- function(DOY, ALAT, TTIMD){
      
      # Compatability issue
      KHRS <- 24
      HHRS <- 12
      
      # Private functions
      ECCENT <- function(T) 0.01675104 - (4.08E-5 + 1.26E-7*T)*T
      
      ANOM <- function(T,D){
        
        anom <- -1.52417 + (1.5E-4 + 3.0E-6*T)*T^2
        anom <- anom + 0.98560*D
        if(anom > 360)
          anom <- anom - 360.0*floor(anom/360.0)
        
        return(anom * pi/180)
      }
      
      EPSIL <- function(T)(23.452294- (1.30125E-2+ (1.64E-6-5.03E-7*T)*T)*T)*pi/180
      
      OMEGA <- function(T,D)(281.22083+ (4.53E-4+3.0E-6*T)*T*T+4.70684E-5*D)*pi/180
      
      ALUNAR <- function(T,D) (259.1833+ (2.078E-3+2.0E-6*T)*T*T-0.05295*D)*pi/180
      # End private functions.
      
      T <- DOY/36525
      
      # COMPUTE SOLAR ORBIT
      ECC <- ECCENT(T)
      
      RM <- ANOM(T,DOY)
      E <- RM
      for(IM in 1:3){
        E <- E + (RM- (E-ECC*sin(E)))/ (1-ECC*cos(E))
      }
      
      V <- 2.0*atan(sqrt((1+ECC)/ (1-ECC))*tan(0.5*E))
      
      if(V < 0) V <- V + 2*pi
      R <- 1 - ECC*cos(E)
      
      EPS <- EPSIL(T)
      OMEG <- OMEGA(T,DOY)
      
      # COMPUTE NUTATION TERMS
      LUNLON <- ALUNAR(T,DOY)
      NUTOBL <- (2.5583E-3+2.5E-7*T)*cos(LUNLON)*pi/180
      EPS <- EPS + NUTOBL
      NUTLON <- - (4.7872E-3+4.7222E-6*T)*sin(LUNLON)*pi/180
      
      # COMPUTE SOLAR DECLINATION
      DEC <- asin(sin(EPS)*sin(V+OMEG))
      
      # COMPUTE EQN OF TIME
      MLON <- OMEG + RM
      if(MLON < 0)MLON <- MLON + 2*pi
      if(MLON > 2*pi)MLON <- MLON - 2*pi*floor(MLON/(2*pi))    
      Y <- (tan(EPS/2))^2
      Y <- (1-Y)/ (1+Y)
      
      SL <- OMEG + NUTLON + V
      if(SL < 0) SL <- SL + 2*pi
      if(SL > 2*pi)SL <- SL - 2*pi*floor(SL/(2*pi))
      AO <- atan(Y*tan(SL))
      
      EQNTIM <- AO - MLON
      EQNTIM <- EQNTIM - pi*floor(EQNTIM/pi)
      if(abs(EQNTIM) > 0.9*pi) EQNTIM <- EQNTIM - pi*EQNTIM/abs(EQNTIM)
      AO <- EQNTIM + MLON
      if(AO > 2*pi) AO <- AO - 2*pi*floor(AO/(2*pi))
      # DAY LENGTH
      MUM <- cos(ALAT-DEC)
      MUN <- -cos(ALAT+DEC)
      MUA <- 0.0
      REFAC <- 0.0
      UMN <- -MUM*MUN
      if(UMN > 0) REFAC = 0.05556/sqrt(UMN)
      if(MUN > MUA) MUA <- MUN    
      if(MUM > MUA){
          FRACSU <- sqrt((MUA-MUN)/ (MUM-MUA))
          FRACSU <- 1.0 - 2.0*atan(FRACSU)/pi
          SUNSET <- HHRS*FRACSU
          SUNRIS <- SUNSET
          SUNSET <- SUNSET + REFAC + EQNTIM*HHRS/pi
          SUNRIS <- SUNRIS + REFAC - EQNTIM*HHRS/pi
          SUNSET <- SUNSET + HHRS + TTIMD
          SUNRIS <- HHRS - SUNRIS + TTIMD
          EQNTIM <- EQNTIM*HHRS/pi
          DAYL <- SUNSET - SUNRIS
      }
    
        
      return(list(DEC=DEC,EQNTIM=EQNTIM,DAYL=DAYL,SUNSET=SUNSET))
      
    }
      
  
  # continue zenaz.
    DATE <- as.Date(ISOdate(year,month,day))
    DJUL <- as.vector(DATE - as.Date("1900-1-1") + 1)
  
	k <- pi/180
	
  # latitude in radians
  ALAT <- lat * k
	
	if(long < 0){
	    long <- 360.0 - long
        tzlong <-  360.0 - tzlong
	}
		
	ALONG <- long * k
  tzlong <- tzlong * k  
	
	if(!LAT){
		TTIMD <- (24/ (2*pi))*(ALONG - tzlong)
	} else {
		TTIMD <- 0
	}

	# Maestra evaluates solar position mid-timestep (see zenaz subroutine).
	if(all(is.na(timeofday))){
		HOURS <- seq(from=24/KHRS/2, by=24/KHRS, length=KHRS)
	} else {
		HOURS <- timeofday
		KHRS <- length(timeofday)
	}
	
	SUNcall <- SUN(DJUL, ALAT, TTIMD)
  
  DEC <- SUNcall$DEC  
  EQNTIM <- SUNcall$EQNTIM
  DAYL <- SUNcall$DAYL
  SUNSET <- SUNcall$SUNSET
  
  # To match old definition, with 24 hours in a day.
  solarnoon <- 12
  
  ZEN <- c()
  AZ <- c()  

  SolarTime <- HOURS - TTIMD - EQNTIM    
      
  # get solar zenith and azimuth angles.
  for(i in 1:KHRS){

      # hour angle
      HI <- (pi/180) * (12 - SolarTime[i])*15

      # zenith angle 
      ZEN[i] <- acos(sin(ALAT)*sin(DEC) + 
        cos(ALAT)*cos(DEC)*cos(HI))
      
      # Cosine of the azimuth angle (Iqbal)
      H <- pi/2 - ZEN[i]
      COSAZ <- (sin(H)*sin(ALAT)- sin(DEC)) /
                   (cos(H)*cos(ALAT))
      if (COSAZ >  1.0)
        COSAZ <- 1.0
      if (COSAZ < -1.0) 
        COSAZ <- -1.0
      
      AZ[i] <- acos(COSAZ)
      
      if(SolarTime[i] > 12)AZ[i] <- -AZ[i]
  }
      
     
  dfr <- data.frame(zen=ZEN, az= pi - AZ)
	dfr[dfr$zen > pi/2,] <- NA
	dfr <- dfr / k
      
	# Solar altitude.
	dfr$alt <- 90 - dfr$zen
	
  return(list(hour=HOURS, altitude=dfr$alt, azimuth=dfr$az, 
	daylength=DAYL, sunset = SUNSET, zenrad=k*dfr$zen, LAT=SolarTime))

}

  
