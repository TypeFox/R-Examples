#'Generate a weather object
#'
#'@description To run Yplant, a weather object needs to be constructed, that contains solar
#'position data, radiation, air temperature, and so on. This function generates
#'a daily diurnal weather dataset using a fairly standard weather generator, or
#'constructs the weather object with user-specified data. See Details.
#'
#'A built-in weather generator simulates the following variables: \describe{
#'\item{altitude,azimuth}{Position of the sun (degrees).}
#'\item{PAR}{Photosynthetically active radiation (mu mol m-2 s-1).}
#'\item{fbeam}{Fraction direct beam of PAR (-).} \item{Tair}{Air temperature
#'(deg C).} \item{VPD}{Vapor pressure deficit.} } The following two variables
#'are user input, and have no within-day variation: \describe{
#'\item{Ca}{Atmospheric CO2 concentration (ppm). (Default = 390ppm).}
#'\item{Patm}{Atmospheric pressure (kPa). (Default = 1.01kPa).} }
#'
#'If you are curious about the algorithm, please check the code (type
#'\code{setMet}).
#'
#'To generate a weather dataset, simply use this command: \preformatted{
#'aprilday <- setMet(richmond, nsteps=12, Tmin=9, Tmax=25, month=6, day=21) }
#'Where \code{richmond} is a Yplant location object, generated with
#'\code{\link{setLocation}}.
#'
#'The weather object can be plotted: the following command produces a simple
#'built-in graph of PAR, Tair, VPD and fbeam: \preformatted{ plot(aprilday) }
#'
#'Alternatively, the user can input a dataframe (or CSV file) that contains the
#'weather variables (or a subset of them). For example, \preformatted{ mymet <-
#'data.frame(Tair=20, PAR0=seq(5,1000,length=10), fbeam=0, Ca=400) } The names
#'of the variables need to be *exactly* as described above (and are
#'case-sensitive!).
#'
#'If solar altitude and azimuth are not provided, they will be calculated from
#'the location object. In the case that \code{fbeam = 0}, though, the solar
#'position has no effect and is ignored, and not calculated.
#'
#'@param location A Yplant location object (class 'yplocation', see
#'\code{\link{setLocation}}).
#'@param metdat Optionally, a dataframe (or name of CSV file) with standard
#'weather variables.
#'@param year Optional (slight effects on solar path).
#'@param month 1-12
#'@param day day of month
#'@param nsteps number of steps (will affect number of simulation steps in
#'\code{\link{YplantDay}}.
#'@param PARday Total daily PAR on a horizontal surface (MJ m-2 d-1).
#'@param AtmTrans Atmospheric transmission.
#'@param fbeamday Daily beam fraction.
#'@param fbeammethod If 'Spitters', uses the Spitters algorithm to estimate
#'fbeam by timestep, otherwise it is constant (and given by \code{fbeamday}.
#'@param Tmin,Tmax Daily minimum and maximum temperature (deg C).
#'@param VPDmax Optional. Daily maximum VPD (if not given, estimated from
#'\code{Tmin}).
#'@param maxlag Lag of temperature maximum behind solar maximum (fraction of
#'day).
#'@param Ca Atmospheric CO2 concentration (ppm).
#'@param Patm Atmospheric pressure (kPa).
#'@return An object of class 'ypmet', a list with the following components:
#'\describe{ 
#'\item{dat}{A dataframe with the weather variables (see Details for
#'a description).} 
#'\item{method}{Either 'generated' (weather generator was
#'used), or 'input' when user provided \code{metdata}.} 
#'\item{daylength}{in hours} 
#'\item{sunset,sunrise}{in hours} 
#'\item{location}{A Yplant location object (class 'yplocation', see \code{\link{setLocation}}} }
#'@author Remko Duursma.  Solar path and diffuse partitioning code borrowed
#'from Maestra (thanks to Belinda Medlyn).
#'@seealso \code{\link{setPhy}},\code{\link{setLocation}}
#'@references For fraction diffuse radiation, uses the 'Spitters algorithm':
#'
#'Spitters, C.J.T., Toussaint, H.A.J.M., Goudriaan, J., 1986, Separating the
#'diffuse and direct component of global radiation and its implications for
#'modeling canopy photosynthesis. Part I. Components of incoming radiation, Ag.
#'For. Meteorol., 38:217-229.
#'@keywords misc
#'@export
setMet <- function(location=NULL, 
                    metdat=NULL,
	                  year=2012,   # not so important; can use default.
        					  month=NA, 
        					  day=NA,  
        					  nsteps=10, 
        					  PARday=22,
        					  AtmTrans=0.76,
        					  fbeamday=NA,
        					  fbeammethod=c("spitters","constant"),
        					  Tmin=10,
        					  Tmax=25,
        					  VPDmax=NA,
        					  maxlag=0.1,
        					  Ca=390,
        					  Patm=101
					      ){
	
	
	fbeammethod <- match.arg(fbeammethod)

	if(!is.na(month) && !is.na(day)){
		z <- zenaz(year, month, day, 
			location$lat, location$long, location$tzlong)
		sunrise <- z$sunset - z$daylength
		sunset <- z$sunset
		dh <- z$daylength / nsteps / 2
		
		# Time of day (in hrs local time).
		hours <- seq(sunrise+dh, sunset-dh, length=nsteps)

		# Zenith, azimuth of sun at many times; then interpolate.
		# Have to take this long route : some problem with fewer timesteps!
		sunposall <- zenaz(year, month, day, 
			location$lat, location$long, location$tzlong,
			timeofday=seq(sunrise+dh, sunset-dh, length=151))
    
		faz <- approxfun(x=sunposall$hour, y=sunposall$azimuth)
		falt <- approxfun(x=sunposall$hour, y=sunposall$altitude)
	}
	
	# Generate weather data, using a simple built-in routine.
	if(is.null(metdat)){
		if(is.null(location))stop("Need at least a location object, or a 'metdat' dataframe.")
		if(!inherits(location, "yplocation"))stop("location should be generated with setLocation().")
		if(is.na(month) || is.na(day))
			stop("Need month and day, to calculate solar path.")
		
		# Find azimuth and altitude.
		Azimuth <- faz(hours)
		Altitude <- falt(hours)
		Zenrad <- pi/2 - Altitude*pi/180

		# Simple estimate of available PAR (Yplant; Pearcy and Yang 1996).
		# a <- sunpos$altitude * pi/180
		# PARnorm <- solcons * AtmTrans^(1/sin(a))
		# PAR <- PARnorm * sin(a)
		DOY <- as.POSIXlt(ISOdate(year,month,day))$yday + 1
		
		pardfr <- calcparhrly(Zenrad, DOY, PARday, nsteps, 
			AtmTrans, fbeamday, fbeammethod)
		
		# VPD and T diurnal.
		vpdt <- VPDTdiurnal(Tmax=Tmax,Tmin=Tmin,reltime=(hours-sunrise)/z$daylength,
							maxlag=maxlag,VPDmax=VPDmax)
		
		dfr <- data.frame(timeofday=hours, altitude=Altitude, 
			azimuth=Azimuth, PAR=pardfr$PAR0, fbeam=pardfr$fbeam, 
			Tair=vpdt$Tair, VPD=vpdt$VPD, Ca=Ca, Patm=Patm)
	}
	if(!is.null(metdat)){
	
		# Some error checking might be nice here...
		if(!is.data.frame(metdat))metdat <- read.csv(metdat)
		nsteps <- nrow(metdat)
		
    if("timeofday" %in% names(metdat))
      hours <- metdat$timeofday
    else
      hours <- 1:nsteps
	
		# If fbeam not 0 (diffuse only), azimuth and altitude not in metdat,
		# OR timeofday AND a location object, we cannot continue!
		if("fbeam" %in% names(metdat))
			diffonly <- all(metdat$fbeam == 0)  # TRUE if all fbeams in metdat = 0.
		else
			diffonly <- FALSE
		
		# Not a setting that would be recommended, but it is possible!
		if(fbeamday==0 & fbeammethod=="constant"){
			diffonly <- TRUE
			metdat$fbeam <- 0
		}
		azimaltset <- ("azimuth" %in% names(metdat)) && ("altitude" %in% names(metdat))
		if(!diffonly & !azimaltset & is.null(location))
			stop("** - Either need altitude&azimuth in metdat, or provide a location object!.")
		
		if(!azimaltset & !diffonly & !("timeofday" %in% names(metdat)))
			stop("Need solar altitude & azimuth in metdat, or provide 'timeofday' to calculate them.")
		if(!azimaltset & !diffonly & (is.na(month) & is.na(day)))
			stop("Need month and day set to calculate solar path.")
		
		if(azimaltset || diffonly){
			# For consistency with print.ypmet.
			z <- list()
			z$daylength <- "not calculated"
			z$sunset <- "not calculated"
			sunrise <- "not calculated"
			
			if(diffonly){
				metdat$azimuth <- 0
				metdat$altitude <- 0
			}
			
			#
			if(!("PAR0" %in% names(metdat)))
				stop("You provided altitude&azimuth, but not PAR0. No current method to estimate PAR0.")
			
		} else {
			# z <- zenaz(year, month, day, 
			# location$lat, location$long, location$tzlong)
			# sunrise <- z$sunset - z$daylength
			# sunset <- z$sunset
			# dh <- z$daylength / nsteps / 2
			
			# Time of day (in hrs local time).
			hours <- seq(sunrise+dh, sunset-dh, length=nsteps)

			# # Zenith, azimuth of sun at those times (list).
			# sunpos <- zenaz(year, month, day, 
				# location$lat, location$long, location$tzlong,
				# timeofday=hours)
			metdat$altitude <- falt(hours) #sunpos$altitude
			metdat$azimuth <- faz(hours) #sunpos$azimuth
			Zenrad <- pi/2 - Altitude*pi/180
			
			# 
			if(!("PAR0" %in% names(metdat)) |  !("fbeam" %in% names(metdat))){
				
				DOY <- as.POSIXlt(ISOdate(year,month,day))$yday + 1
				pardfr <- calcparhrly(Zenrad, DOY, PARday, nsteps, 
				AtmTrans, fbeamday, fbeammethod)
				if(!("PAR0" %in% names(metdat))){
					message("Calculating incident PAR (because PAR0 not given).")
					metdat$PAR0 <- pardfr$PAR0
				}
				if(!("fbeam" %in% names(metdat))){
					message("Calculating beam fraction (because fbeam not given)")
					metdat$fbeam <- pardfr$fbeam
				}	
			}
		}
		
		# VPD and T diurnal. Do NOT calculate if not given (too much trouble).
		if(!("VPD" %in% names(metdat))){
			metdat$VPD <- 1.5
			message("Setting VPD=1.5 kPa, because VPD not in met data.")
		}
		if(!("Tair" %in% names(metdat))){
			message("Setting Tair=25 degC, because Tair not in met data.")
			metdat$Tair <- 25
		}
		if(!("Ca" %in% names(metdat)))
			metdat$Ca <- Ca
		if(!("Patm" %in% names(metdat)))
			metdat$Patm <- Ca

		dfr <- data.frame(timeofday=hours, altitude=metdat$altitude, 
			azimuth=metdat$azimuth, PAR=metdat$PAR0, fbeam=metdat$fbeam, 
			Tair=metdat$Tair, VPD=metdat$VPD, Ca=metdat$Ca, Patm=metdat$Patm)

	}
	
l <- list()
l$dat <- dfr
l$method <- ifelse(is.null(metdat),"generated","input")
l$daylength <- z$daylength
l$sunset <- z$sunset
l$sunrise <- sunrise
l$location <- location
l$year <- year
l$month <- month
l$day <- day
class(l) <- "ypmet"	


return(l)	
}

