astrocalc4r <- function(day,month,year,hour,timezone,lat,lon,withinput=FALSE,seaorland="maritime",acknowledgment=FALSE)
  {
  #version 2.2 AstroCalcPureRv
  if (acknowledgment){
    cat("\n","---------------------------------------------------------")
    cat("\n","                AstroCalcPureR Version 2.2")
    cat("\n","Documentation: Jacobson L, Seaver A, Tang J. 2011. AstroCalc4R:")
    cat("\n","software to calculate solar zenith angle; time at sunrise,")
    cat("\n","local noon and sunset; and photosynthetically available")
    cat("\n","radiation based on date, time and location. US Dept Commer,")
    cat("\n","Northeast Fish Sci Cent Ref Doc. 11-14; 10 p. Available from:")
    cat("\n","National Marine Fisheries Service, 166 Water Street, ")
    cat("\n","Woods Hole, MA 02543-1026, or online at")
    cat("\n","http://nefsc.noaa.gov/publications/")
    cat("\n \n","Code available gratis at http://www.nefsc.noaa.gov/AstroCalc4R/")
    cat("\n \n","Contact: Larry.Jacobson@noaa.gov")
	cat("\n\n","Useage:")
	cat("\n","    AstroCalcPureR(day,month,year,hour,timezone,")
	cat("\n",  "                   lat,lon,withinput=F,")
	cat("\n",  "                   seaorland='maritime',")
	cat("\n",  "                   acknowledgment=TRUE)")
	cat("\n\n","HINT: set acknowledgment=FALSE to avoid this message")
    cat("\n","---------------------------------------------------------","\n")
  }  
         
  options(digits=9)
  deg2rad <- pi/180

####
# Check for null or missing data vectors
  null.c<-function(x) return(sum(is.null(x))) 
	if (sum(null.c(day),null.c(month),null.c(year),null.c(hour),null.c(timezone),null.c(lat),null.c(lon)) > 0) 
		stop("\n Null or missing required data vector for day, month, year, timezone, lat or lon \n")
		  
####
# Be sure all arguments have the same length
	if ((length(day)!=length(month))|(length(month)!=length(year))|
		(length(year)!=length(hour))|(length(hour)!=length(timezone))|
		(length(timezone)!=length(lat))|(length(lat)!=length(lon)))
	stop("\n Input vectors are not the same length \n")
	
# Save n records	
	times<-length(day)
	
####
# Check for NA values
	na.c <- function(x) return(sum(is.na(x))) #simple function to simplify 
	if (sum(na.c(day),na.c(month),na.c(year),na.c(hour),na.c(timezone),na.c(lat),na.c(lon)) > 0) 
		stop("\n NA values in input data \n")

####		
# Check for obviously incorrect data 

# Month and date (considering leap years)
	logic1 <- year < 0
	if (sum(logic1)>0) stop(cat("\n Error in year for record numbers:",(1:times)[logic1]," \n\n"))

	is.leap <- function(x) return ((((x %% 4==0)&(x %% 100 !=0)))|(x %% 400 ==0))
	date.list <- c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
	
	logic1 <- abs(month-6) > 6
	if (sum(logic1) > 0) stop(cat("\n Error in month for record numbers:",(1:times)[logic1]," \n\n"))
	
	logic1 <- day > (date.list[month] + is.leap(year)*(month==2))
	logic2 <- day <= 0
	if ((sum(logic1) > 0) |(sum(logic2) >0)) stop(cat("\n Incorrect month-day-year combination for record numbers: ",(1:times)[logic1|logic2]," \n\n"))
	
	logic1 <- abs(hour-12) > 12
	if (sum(logic1) > 0) stop(cat("\n Error in hour for record numbers:",(1:times)[logic1]," \n\n"))
	
	logic1 <- abs(timezone)>12
	if (sum(logic1)>0) stop(cat("\n Error in time zone for record numbers:",(1:times)[logic1]," \n\n"))

	logic1 <- abs(lon) > 180
	if (sum(logic1)>0) stop(cat("\n Error in longitude for record numbers:",(1:times)[logic1]," \n\n"))
	
	logic1 <- abs(lat) > 90
	if (sum(logic1)>0) stop(cat("\n Error in latitude for record numbers:",(1:times)[logic1]," \n\n"))
	
  
#Calculate Julian Day Starting at 4712 BC
#Citation: "Astronomical Algorithms"  P.61
	JulianDay <- function(xday,xmonth,xyear)
		{
		mm <- xmonth

		xmonth[mm <=2] <- xmonth[mm <=2] +12
		xyear[mm <=2] <- xyear[mm<=2] -1
		xa <- floor(xyear/100)
		xb <- 2 - xa + floor(xa/4)
		jd <- floor(365.25*(xyear + 4716)) + floor(30.6001 * (xmonth + 1)) + xday + xb - 1524.5
		return(jd)
		}

	daymonth <- function(mth,yr)
		{
		day[is.leap(yr)] <- c(31,29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)[mth[is.leap(yr)]]
		day[!is.leap(yr)] <- c(31,28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)[mth[!is.leap(yr)]]
		return(day)
		}
 
#**  Calculate Photosynthetically Available Radiation (PAR - Watts / M2)
#**
#**  Basis: Fouin, et al "Journal of Geophysical Research, Vol 94, No. C7 pp. 9731-9742 
#**  July 15, 1989  
#**
#**  Equation 6
#**
#**  equation 1 in frouin et al. 1989 paper should be (d0/d)2 not as it is in the paper (d/d0)2  
#**  in their equations 1 and 6
#**  email from r.frouin to j.o'reilly, oct 28,2003:
#**    "there is a typo in the equation, it's actually (do/d)2. obviously when the actual distance (d) is larger than the average distance
#**    "(do), the correction factor should be smaller (less irradiance since the sun is farther away)."
#**
#*/
	parcalc <- function(zenith,setting=seaorland)
		{
#NOTE: this function is now vectorized Jiashen Tang 04/05/2013 

# All calculations are for PAR at 300-700 nm (the most popular definition)
#  but 350-700 and 250-4000 nm can be implemented using the parameters 
#  in Table 2 of Frouin, R., Lingner, D., Gautier, C., Baker, K. 
#  and Smith, R.  1989.  A simple analytical formula to compute total 
#  and photosynthetically available solar irradiance at the ocean surface
#  under clear skies.  J. Geophys. Res. 94: 9731-9742. 

#Maritime 400-700nm

		I0 <- 531.2 #Solar Irradiance
		V <- 23
		uv <- 1.4
		u0 <- .34
		r <- .05
		d <- 1.0
		
		if (!setting %in% c("maritime","continental")) 
			stop("setting value is neither 'maritime' nor 'continental'!")
			
		if (setting=="maritime")
		# in case aerosol type is maritime
			{
			a <- .068
			b <- .379
			a1 <- .117
			b1 <- .493
			av <- .002
			bv <- .87
			a0 <- .052
			b0 <- .99
			} else if (setting == "continental")
		# or aerosol type continental
			{
			a <- .078
			b <- .882
			a1 <- .123
			b1 <- .594
			av <- .002
			bv <- .87
			a0 <- .052
			b0 <- .99
			}
		
			zrad <- zenith * deg2rad
			x1 <- uv / cos(zrad)
			xx <- exp(-av * x1^bv)
			x2 <- u0 / cos(zrad)
			xxx <- exp(-a0 * x2^b0)
			xa <- a +b /V
			xb <- d - r *(a1 + b1/V)
			par <- I0 * cos(zrad) * exp(-xa/cos(zrad))/xb*xx*xxx
			par[zenith >89.9999 ] <- 0
			return(par)
		
		}
	
	output <- as.data.frame(matrix(nrow=0,ncol=9))
	names(output) <- c("noon","sunrise","sunset","azimuth","zenith","eqtime","declin","daylight","PAR" )
	
	
		hourtemp <- hour-timezone

		hour <- ifelse(hourtemp >24,hourtemp-24,hourtemp)
		change_day <- !(hour==hourtemp)

		dm <- daymonth(month,year)
		daytemp <- day

		daytemp[change_day] <- ifelse((day[change_day] < dm[change_day]), day[change_day]+1,1)

		change_month <- abs(day-daytemp) >1

		monthtemp <- month
		monthtemp[change_month] <- ifelse(month[change_month] < 12,month[change_month]+1,1)

		change_year <- abs(month-monthtemp) >1
		yeartemp <- year
		yeartemp[change_year] <- year[change_year]+1


		xy <- yeartemp
		xm <- monthtemp
		xd <- daytemp + hourtemp/24
		# Julian day[i]
		jd <- JulianDay(xd,xm,xy)*100/100
		#Julian Century; "Astro algo" Eq.25.1
		jc <- (jd-2451545)/36525

		#Geometric Mean Longitude of the Sun; Eq. 25.2
		xx <- 280.46646 + 36000.76983*jc + 0.0003032*jc^2
		gmls <- xx %% 360  #L0

		# Mean Anomaly of the Sun; Eq. 25.3
		xx <- 357.52911 + 35999.05029* jc - 0.0001537*jc^2
		gmas <- xx %% 360 #M

		#Eccentricity of the Earth's orbit; Eq. 25.4
		eeo <- .016708634 - .000042037*jc - .0000001267*jc^2

		#Sun's Equation of the Center; p.164
		scx <- (1.914602 - .004817*jc - .000014*jc^2)*sin(gmas*deg2rad) +
		  (.019993 - .000101*jc)*sin(2*gmas*deg2rad)+
		  .000289*sin(3*gmas*deg2rad)
		  
		#Sun's True Longitude & Anomaly; p.164
		Stl <- gmls + scx
		Sta <- gmas + scx
		#the Sun's Radius Vector; Eq. 25.5
		srv <- 1.000001018*(1-eeo^2)/(1 + eeo*cos(Sta*deg2rad)) #R

		#the Sun's Apparent Longitude; p.164
		omega <- 125.04 - 1934.136 *jc
		lambda <- Stl - .00569 - .00478*sin(omega*deg2rad)

		# Mean Obliquity of the Ecliptic; Eq. 22.2
		epsilon <- (23 + 26/60 + 21.448/60^2) - (46.815/60^2)*jc - (.00059/60^2)*jc^2 + (.001813/60^2)*jc^3

		# Obliquity Correction; Eq. 25.8
		oblx <- .00256 * cos(omega*deg2rad)
		epsilon <- epsilon + oblx

		# the Sun's Right Ascension; Eq. 25.6
		#Note: apparent position
		alpha <- atan2(cos(epsilon*deg2rad)*sin(lambda*deg2rad), cos(lambda*deg2rad))/deg2rad

		# the Sun's Declination; Eq. 25.7
		declin <- asin(sin(epsilon*deg2rad)*sin(lambda*deg2rad))/deg2rad

		# Equation of Time; Eq. 28.3
		y <- tan(epsilon*deg2rad/2)^2
		eqtime <- (y*sin(2*gmls*deg2rad) - 2*eeo* sin(gmas*deg2rad) + 4*eeo*y*sin(gmas*deg2rad)*cos(2*gmls*deg2rad)
		      - y^2*sin(4*gmls*deg2rad)/2 - 5/4*eeo^2*sin(2*gmas*deg2rad))/deg2rad*4
		      
		#hour[i] Angle; Eq. 15.1
		h0 <-  -.8333*deg2rad
		phi <- lat*deg2rad
		
		hangle  <- acos((sin(h0)-sin(declin*deg2rad)*sin(phi))/cos(declin*deg2rad)/cos(phi))/deg2rad
			   #sin(h0*deg2rad)/cos(declin*deg2rad)/cos(lat[i]*deg2rad)-tan(declin*deg2rad)* tan(lat[i]*deg2rad)
			   
		#Solar Noon(LST) in minutes
		#Each 15 Degrees of Longitude = 1 hour[i]
		#Each Time zone = 1 hour[i]
		#1440 Minutes in day[i]
		noon <- (720 - 4*lon + timezone*60 - eqtime) /1440

		# Sunrise & Sunset
		sunrise <- (noon *1440 -hangle *4)/1440*24
		sunset <- (noon *1440 + hangle *4)/1440*24
		noon <- noon * 24

		#Length of day
		daylight <- hangle * 8

		# True Solar Time (minutes)
		tst <- (hourtemp*60 + eqtime + 4*lon) %%1440

		# True Solar Angle (degrees)
		#if (tst <0) tsa <- tst /4 +180 else tsa <- tst/4 -180 vectorized
		tsa <- ifelse(tst<0,tst /4 +180,tst/4 -180)
		

		#Zenith; Eq. 13.6
		zenith <- 90 - asin(sin(lat*deg2rad) * sin(declin*deg2rad)+cos(lat*deg2rad)*cos(declin*deg2rad)*cos(tsa*deg2rad))/deg2rad

		#Azimuth: p.94
		azimuth <- acos((sin(lat*deg2rad)*sin((90-zenith)*deg2rad)-sin(declin*deg2rad))/cos(lat*deg2rad)/cos((90-zenith)*deg2rad))/deg2rad+180

		#if (tsa >0) azimuth <- azimuth %%360 else azimuth <- 360 - azimuth %%360
		azimuth <-ifelse(tsa >0,azimuth %%360,360 - azimuth %%360)
		
		
		daylight <- daylight /60
		
		PAR <- parcalc(zenith)
		
		# deal with the Polar situ
		if (any(is.nan(sunrise))) 
			{
			
			message(paste("Warning: Polar day/night (daylength 0 or 24 hrs) at record(s):",(1:times)[is.nan(sunrise)],
							"\n Check input data (i.e. latitude)?")) 
		
			daylight <- ifelse(PAR>0,24,0)
			
			}
		output <- rbind(output,data.frame(noon=noon,sunrise=sunrise,sunset=sunset,azimuth=azimuth,
					zenith=zenith,eqtime=eqtime,declin=declin,
					daylight=daylight,PAR=PAR))
					
#		print(paste("output in round",i,"is:"))
#		print(output)			
# No more for loop so the ending bracket		}
	

	if (withinput) return(cbind(data.frame(tzone=timezone,day=day,month=month,year=year,hhour=hour,xlat=lat,xlon=lon),output)) else return(output)

	}
