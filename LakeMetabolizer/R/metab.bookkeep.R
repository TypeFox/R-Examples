#'@title Metabolism model based on simple day/night summation NEP-interpreted changes in DO.
#'@description This model is a simple model based on the assumption that movements in DO during 
#'the day are due to NEP and gas exchange. Respiration is estimated from night-time decreases. 
#'GPP is calculated from the algebraic manipulation of NEP and R. Based on Cole et al 2000.
#'@param do.obs Vector of dissovled oxygen concentration observations, mg L^-1
#'@param do.sat Vector of dissolved oxygen saturation values based on water temperature. Calculate using \link{o2.at.sat}
#'@param k.gas Vector of kGAS values calculated from any of the gas flux models 
#'(e.g., \link{k.cole}) and converted to kGAS using \link{k600.2.kGAS}
#'@param z.mix Vector of mixed-layer depths in meters. To calculate, see \link{ts.meta.depths}
#'@param irr Integer vector of 1's (daytime) and 0's (nighttime), or numeric vector of irradiance that will be converted to boolean 1's and 0's if "datetime" is passed via \code{...}
#'@param ... additional arguments to be passed, particularly \code{POSIXct} class "datetime"
#'@return
#'A data.frame with columns corresponding to components of metabolism 
#'\describe{
	#'\item{GPP}{numeric estimate of Gross Primary Production, \eqn{mg O_2 L^{-1} d^{-1}}{mg O2 / L / d}}
	#'\item{R}{numeric estimate of Respiration, \eqn{mg O_2 L^{-1} d^{-1}}{mg O2 / L / d}}
	#'\item{NEP}{numeric estimate of Net Ecosystem production, \eqn{mg O_2 L^{-1} d^{-1}}{mg O2 / L / d}}
#'}
#'@author
#'R. Iestyn Woolway, Hilary Dugan, Luke A Winslow, Ryan Batt, Jordan S Read, GLEON fellows
#'@references
#'Cole, Jonathan J., Michael L. Pace, Stephen R. Carpenter, and James F. Kitchell. 2000. 
#'\emph{Persistence of Net Heterotrophy in Lakes during Nutrient Addition and Food Web Manipulations}. 
#'Limnology and Oceanography 45 (8): 1718-1730. doi:10.4319/lo.2000.45.8.1718.
#'@seealso
#'\link{metab.bayesian}, \link{metab.mle}, \link{metab.kalman}
#'@examples
#'library(rLakeAnalyzer)
#'Sys.setenv(TZ='GMT')
#'
#'doobs = load.ts(system.file('extdata', 
#'                            'sparkling.doobs', package="LakeMetabolizer"))
#'wtr = load.ts(system.file('extdata', 
#'                          'sparkling.wtr', package="LakeMetabolizer"))
#'wnd = load.ts(system.file('extdata', 
#'                          'sparkling.wnd', package="LakeMetabolizer"))
#'
#'#Subset a day
#'mod.date = as.POSIXct('2009-07-08', 'GMT')
#'doobs = doobs[trunc(doobs$datetime, 'day') == mod.date, ]
#'wtr = wtr[trunc(wtr$datetime, 'day') == mod.date, ]
#'wnd = wnd[trunc(wnd$datetime, 'day') == mod.date, ]
#'
#'k.gas = k600.2.kGAS.base(k.cole.base(wnd[,2]), wtr[,3], 'O2')
#'do.sat = o2.at.sat.base(wtr[,3], altitude=300)
#'
#'# Must supply 1 for daytime timesteps and 0 for nighttime timesteps
#'irr = as.integer(is.day(doobs[,1], 45))
#'
#'metab.bookkeep(doobs[,2], do.sat, k.gas, z.mix=1, irr, datetime=doobs$datetime)
#'@export
metab.bookkeep <- function(do.obs, do.sat, k.gas, z.mix, irr, ...){

  complete.inputs(do.obs=do.obs, do.sat=do.sat, k.gas=k.gas, 
                  z.mix=z.mix, irr=irr, error=TRUE)
  
	nobs <- length(do.obs)  

	mb.args <- list(...)
  
	if(any(z.mix <= 0)){
	  stop("z.mix must be greater than zero.")
	}
  
	if("datetime"%in%names(mb.args)){ # check to see if datetime is in the ... args
		datetime <- mb.args$datetime # extract datetime
		freq <- calc.freq(datetime) # calculate sampling frequency from datetime
		if(nobs!=freq){ # nobs and freq should agree, if they don't issue a warning
			bad.date <- format.Date(datetime[1], format="%Y-%m-%d")
			warning("number of observations on ", bad.date, " (", nobs, ") ", "does not equal estimated sampling frequency", " (", freq, ")", sep="")
		}
	}else{ # if datetime is *not* in the ... args
		warning("datetime not found, inferring sampling frequency from # of observations") # issue a warning (note checks in addNAs)
		# NOTE: because of the checks in addNA's, it is unlikely a user would receive this warning via metab()
		# warning will only be seen through direct use of metab.bookkeep when datettime is not supplied
		freq <- nobs
	}


	if(all(c("datetime", "lake.lat")%in%names(mb.args))){ # check to see if datetime and lake.late are in mb.args
		irr <- as.integer(is.day(datetimes=datetime, lat=mb.args$lake.lat)) # calculate 1's and 0's for irr from datetime and lake.lat
		dayI <- irr == 1L # TRUE when the lights are on
		nightI <- irr == 0L	# TRUE when the lights are off
	}else{
		if(!all(irr==1L | irr==0L)){ # datetime/lake.lat are not found, and if all elements of irr are not integer 0/1
			stop("either supply datetime & lake.lat arguments, or supply irr as integer vector of 1's and 0's") # then that's an error
		}
		# but if datetime/lake.lat are not found in mb.args, yet irr is all integer 0/1, then ...
		dayI <- irr == 1L 
		nightI <- irr == 0L
	}

	delta.do <- diff(do.obs)
	miss.delta <- sum(is.na(delta.do)) # number of NA's
	if(miss.delta != 0){ # the number of NA's should be 0, otherwise issue a warning
		warning(paste(miss.delta, " missing values (", miss.delta/length(delta.do), "%) in diff(do.obs)", sep=""))
		# Note: it should be hard to get this warning. Should be impossible to get this warning through metab due to addNA's
	}

	#gas flux out is negative
	#normalized to z.mix, del_concentration/timestep (e.g., mg/L/10min)
	gas.flux <- (do.sat - do.obs) * (k.gas/freq) / z.mix 

	#remove the component of delta.do that is due to gas flux
	delta.do.metab <- delta.do + gas.flux[1:(length(gas.flux)-1)]

	#normalize units to per-day
	# delta.do.meta.daily <- delta.do.metab * (60*60*24)/as.numeric(delta.times, 'secs')

	nep.day <- delta.do.metab[dayI]
	nep.night <- delta.do.metab[nightI]


	R <- mean(nep.night, na.rm=TRUE) * freq # should be negative
	NEP <- mean(delta.do.metab, na.rm=TRUE) * freq # can be positive or negative
	GPP <- mean(nep.day, na.rm=TRUE) * sum(dayI) - R # should be positive

	metab <- data.frame("GPP"=GPP, "R"=R, "NEP"=NEP)
	return(metab)

}

