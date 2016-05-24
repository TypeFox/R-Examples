#' @title Calculates the equilibrium saturation concentration of oxygen in water
#'   at the supplied conditions
#' @description Used to calculate the equilibrium concentration of oxygen in 
#'   water. The equilibration concentration of oxygen in water varies with both 
#'   temperature, salinity, and the partial pressure of oxygen in contact with 
#'   the water (calculated from supplied elevation or barometric pressure).
#' @details DO solubility is converted from mL/L to mg/L by multiplying by 
#'   1.42905, per USGS memo 2011.03. Corrections for vapor pressure are made 
#'   according to barometric pressure as in Equations 2&3 of USGS memos 81.11 
#'   and 81.15. When barometric pressure is not supplied, it is estimated from 
#'   altitude by the barometric formula as in Colt (2012).
#' @name o2.at.sat
#' @param ts.data Object of class data.frame with two named columns 
#'   \dQuote{datetime} and \dQuote{wtr} (water temp in deg C).
#' @param temp a numeric vector of water temperature in degrees Celsius.
#' @param baro barometric pressure in millibars.
#' @param altitude a numeric value indicating the elevation above mean sea level
#'   in meters. Defaults to mean sea level. An alternative to supplying 
#'   barometric pressure.
#' @param salinity a numeric vector of salinity in PSU. Defaults to zero. Length
#'   must be one or equal to length of temperature.
#' @param model the empirical model to be used. \code{"garcia-benson"}, 
#'   \code{"garcia"}, \code{"weiss"} and \code{"benson"} are the available 
#'   options. \code{"garcia-benson"} is our current recommendation. The models 
#'   correspond to the like-named references described below, where both 
#'   \code{"garcia"} and \code{"garcia-benson"} are from Garcia & Gordon (1992).
#'   
#' @return The equilibration concentration at the supplied conditions in mg/L of
#'   oxygen.
#' @author Luke A Winslow
#' @references
#' 
#' Colt, John. \emph{1 - Solubility of Atmospheric Gases in Freshwater.} In 
#' Computation of Dissolved Gas Concentration in Water as Functions of 
#' Temperature, Salinity and Pressure (Second Edition), edited by John Colt, 
#' 1-71. London: Elsevier, 2012.
#' http://www.sciencedirect.com/science/article/pii/B9780124159167000012.
#' 
#' Garcia, H., and L. Gordon (1992), \emph{Oxygen solubility in seawater: Better
#' fitting equations}, Limnol. Oceanogr., 37(6).
#' 
#' Benson, B. B. & Krause, D. (1984). \emph{The concentration and isotopic 
#' fractionation of oxygen dissolved in freshwater and seawater in equilibrium 
#' with the atmosphere.} Limnology and Oceanography, 29(3), 620-632. 
#' doi:10.4319/lo.1984.29.3.0620
#' 
#' Staehr, Peter A., Darren Bade, Matthew C. Van de Bogert, Gregory R. Koch,
#' Craig Williamson, Paul Hanson, Jonathan J. Cole, and Tim Kratz. \emph{Lake
#' Metabolism and the Diel Oxygen Technique: State of the Science.} Limnology
#' and Oceanography: Methods 8, no. 11 (November 1, 2010): 628-44.
#' doi:10.4319/lom.2010.8.0628
#' 
#' USGS. \emph{New Tables of Dissolved Oxygen Saturation Values.} Quality of 
#' Water Branch, 1981. http://water.usgs.gov/admin/memo/QW/qw81.11.html.
#' 
#' USGS. \emph{New Tables of Dissolved Oxygen Saturation Values; Amendment of 
#' Quality of Water Technical Memorandum No. 81.11.} Quality of Water Branch, 
#' 1981. http://water.usgs.gov/admin/memo/QW/qw81.15.html.
#' 
#' USGS. \emph{Change to Solubility Equations for Oxygen in Water.} Technical 
#' Memorandum 2011.03. USGS Office of Water Quality, 2011.
#' 
#' Weiss, R. (1970). \emph{The solubility of nitrogen, oxygen and argon in water
#' and seawater}. Deep Sea Research and Oceanographic Abstracts, 17(4), 721-735.
#' doi:10.1016/0011-7471(70)90037-9
#' 
#' @seealso \link{water.density}, \link{o2.at.sat.base}
#' @keywords math, methods
#' @examples
#' temp.range = 1:25
#' sal.range = 1:25
#' 
#' par(mfrow=c(1,2))
#' plot(temp.range, o2.at.sat.base(temp.range), xlab='Temperature (C)', 
#' ylab='Oxygen Saturation (mg/L)')
#' plot(o2.at.sat.base(rep(20,25), salinity=sal.range), xlab='Salinity (PSU)', ylab='')
#' 
#' @export
o2.at.sat <- function(ts.data, baro, altitude=0, salinity=0, model='garcia-benson'){
	if(ncol(ts.data) > 2){
		stop('Temp can only have two columns, "datetime" and "wtr"')
	}

	dosat <- o2.at.sat.base(ts.data$wtr, baro, altitude, salinity, model)

	return(data.frame(datetime=ts.data$datetime, do.sat=dosat))  
}


#' @rdname o2.at.sat
#' @export
o2.at.sat.base <- function(temp, baro, altitude=0, salinity=rep(0,length(temp)), model='garcia-benson'){
  
  # Conversion from mL/L (the usual output of the garcia, weiss, etc. equations)
  # to mg/L per USGS memo 2011.03
  mgL.mlL <- 1.42905

  # Correction for air pressure; incorportes effects of altitude & vapor pressure of water
  mmHg.mb <- 0.750061683 # conversion from mm Hg to millibars
  if(missing(baro)){
    mmHg.inHg <- 25.3970886 # conversion from inches Hg to mm Hg
    standard.pressure.sea.level <- 29.92126 # Pb, inches Hg
    standard.temperature.sea.level <- 15 + 273.15 # Tb, 15 C = 288.15 K
    gravitational.acceleration <- 9.80665 # g0, m/s^2
    air.molar.mass <- 0.0289644 # M, molar mass of Earth's air (kg/mol)
    universal.gas.constant <- 8.31447 #8.31432 # R*, N*m/(mol*K)
    
    # estimate pressure by the barometric formula
    baro <- (1/mmHg.mb) * mmHg.inHg * standard.pressure.sea.level * 
	    exp( (-gravitational.acceleration * air.molar.mass * altitude) / (universal.gas.constant * standard.temperature.sea.level) )
	}
  # pressure correction per USGS memos 81.11 and 81.15. calculate u by Antoine equation.
  u <- 10 ^ (8.10765 - 1750.286 / (235 + temp)) # u is vapor pressure of water; water temp is used as an approximation for water & air temp at the air-water boundary
  press.corr <- (baro*mmHg.mb - u) / (760 - u) # pressure correction is ratio of current to standard pressure after correcting for vapor pressure of water. 0.750061683 mmHg/mb
  
  # Estimate O2 at saturation in mL/L by several models
	if(tolower(model) == 'garcia'){

	  Ts <- log((298.15 - temp)/(273.15 + temp))

	  lnC <- 2.00856 + 3.22400 *Ts + 3.99063*Ts^2 + 4.80299*Ts^3 + 9.78188e-1*Ts^4 + 
	    1.71069*Ts^5 - salinity*(6.24097e-3 + 6.93498e-3*Ts + 6.90358e-3*Ts^2 + 4.29155e-3*Ts^3) - 3.1168e-7*salinity^2

	  o2.sat <- exp(lnC)

	} else if(tolower(model) == 'garcia-benson'){
	  
	  Ts <- log((298.15 - temp)/(273.15 + temp))
	  
	  lnC <- 2.00907 + 3.22014*Ts + 4.05010*Ts^2 + 4.94457*Ts^3 + -2.56847e-1*Ts^4 + 
	    3.88767*Ts^5 - salinity*(6.24523e-3 + 7.37614e-3*Ts + 1.03410e-2*Ts^2 + 8.17083e-3*Ts^3) - 4.88682e-7*salinity^2
	  
	  o2.sat <- exp(lnC)
	  
	} else if(tolower(model) == 'weiss'){
		tempk <- temp + 273.15

		lnC <- -173.4292 + 249.6339 * (100 / tempk) + 143.3483 *
		log(tempk / 100) - 21.8492 * (tempk / 100) + 
		salinity * (-0.033096 + 0.014259 * (tempk / 100) - 0.0017000 * (tempk / 100)^2)
          
		o2.sat <- exp(lnC)

	} else if(tolower(model) == 'benson'){
		## TODO: Fix this to include salinity
		if(!all(salinity==0)){
			warning('Benson model does not currently include salinity')
		}

	  o2.sat <- (-0.00006 * (temp)^3) + (0.00725 * (temp)^2) - (0.39571 * (temp)) + 14.59030
	  
	  o2.sat <- o2.sat / mgL.mlL # undo the conversion (below) from ml/L to mg/L; Benson model appears to already predict in mg/L
	  
	} else {
	  stop(paste0('unrecognized model: ', model))
	}
	o2.sat <- o2.sat * mgL.mlL * press.corr

	return(o2.sat)
}
