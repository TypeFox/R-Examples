#'@title Metabolism model based on a ordinary least squares parameter estimation framework.
#'@description This function runs the ordinary least squares metabolism model on the supplied gas concentration and other supporting data. This is a 
#'common approach that allows for the concurrent estimation of metabolism paramters from a timeseries.
#'@param do.obs Vector of dissolved oxygen concentration observations, mg L^-1
#'@param do.sat Vector of dissolved oxygen saturation values based on water temperature. Calculate using \link{o2.at.sat}
#'@param k.gas Vector of kGAS values calculated from any of the gas flux models 
#'(e.g., \link{k.cole}) and converted to kGAS using \link{k600.2.kGAS}
#'@param z.mix Vector of mixed-layer depths in meters. To calculate, see \link{ts.meta.depths}
#'@param irr Vector of photosynthetically active radiation in \eqn{\mu mol\ m^{-2} s^{-1}}{micro mols / m^2 / s}
#'@param wtr Vector of water temperatures in \eqn{^{\circ}C}{degrees C}. Used in scaling respiration with temperature
#'@param ... additional arguments; currently "datetime" is the only recognized argument passed through \code{...}
#'@return
#'A data.frame with columns corresponding to components of metabolism 
#'\describe{
	#'\item{GPP}{numeric estimate of Gross Primary Production, \eqn{mg O_2 L^{-1} d^{-1}}{mg O2 / L / d}}
	#'\item{R}{numeric estimate of Respiration, \eqn{mg O_2 L^{-1} d^{-1}}{mg O2 / L / d}}
	#'\item{NEP}{numeric estimate of Net Ecosystem production, \eqn{mg O_2 L^{-1} d^{-1}}{mg O2 / L / d}}
#'}
#'
#'@importFrom stats lm model.matrix
#'
#'@author Luke A Winslow, Ryan Batt, GLEON Fellows
#'@seealso
#'\link{metab}, \link{metab.bookkeep}, \link{metab.mle}, \link{metab.kalman}, \link{metab.bayesian}, 
#'@examples
#'library(rLakeAnalyzer)
#'doobs = load.ts(system.file('extdata', 
#'                            'sparkling.doobs', package="LakeMetabolizer"))
#'wtr = load.ts(system.file('extdata', 
#'                          'sparkling.wtr', package="LakeMetabolizer"))
#'wnd = load.ts(system.file('extdata', 
#'                          'sparkling.wnd', package="LakeMetabolizer"))
#'irr = load.ts(system.file('extdata', 
#'                          'sparkling.par', package="LakeMetabolizer"))
#'
#'#Subset a day
#'mod.date = as.POSIXct('2009-07-08')
#'doobs = doobs[trunc(doobs$datetime, 'day') == mod.date, ]
#'wtr = wtr[trunc(wtr$datetime, 'day') == mod.date, ]
#'wnd = wnd[trunc(wnd$datetime, 'day') == mod.date, ]
#'irr = irr[trunc(irr$datetime, 'day') == mod.date, ]
#'z.mix = ts.thermo.depth(wtr)
#'
#'k600 = k.cole.base(wnd[,2])
#'k.gas = k600.2.kGAS.base(k600, wtr[,3], 'O2')
#'do.sat = o2.at.sat.base(wtr[,3], altitude=300)
#'
#'metab.ols(doobs[,2], do.sat, k.gas, z.mix[,2], irr[,2], wtr[,3])
#'@export
metab.ols <- function(do.obs, do.sat, k.gas, z.mix, irr, wtr, ...){

  complete.inputs(do.obs=do.obs, do.sat=do.sat, k.gas=k.gas, 
                  z.mix=z.mix, irr=irr, wtr=wtr, error=TRUE)
  
	nobs <- length(do.obs)
	
	mo.args <- list(...)
  
	if(any(z.mix <= 0)){
	  stop("z.mix must be greater than zero.")
	}
	if(any(wtr <= 0)){
		stop("all wtr must be positive.")
	}
  
	if("datetime"%in%names(mo.args)){ # check to see if datetime is in the ... args
		datetime <- mo.args$datetime # extract datetime
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

	do.diff <- diff(do.obs)


	inst_flux <- (k.gas/freq) * (do.sat - do.obs)  # positive is into the lake
	
	# flux <- (inst_flux[1:(n.obs-1)] + inst_flux[-1])/2
	flux <- inst_flux[-nobs]

	
	noflux.do.diff <- do.diff - flux/z.mix[-nobs]
	
	lntemp <- log(wtr)
	mod <- lm(noflux.do.diff ~ irr[-nobs] + lntemp[-nobs] -1) # note that if we decied to use resp=rho*log(Temp), you would do lm(do~irr+lntemp-1) (mind the -1)

	rho <- mod[[1]][2] 
	iota <- mod[[1]][1]
	mod.matrix <- model.matrix(mod)
	gpp <- mean(iota*mod.matrix[,1], na.rm=TRUE) * freq
	resp <- mean(rho*mod.matrix[,2], na.rm=TRUE) * freq
	nep <- gpp + resp
	#other ways to get nep:
	# nep = rho + gpp
	# nep = sum(fitted(mod), na.rm=TRUE) # even if there are NA's in the response variable, they shouldn't be included in fitted() ....
	# sum(noflux.do.diff) # i think this should be the same as sum(fitted(mod)) b/c model residuals sum to 0 ... right?
	# also note that NEP is gpp+rho (rho is negative by this convention, which is consistent w/ Kalman, Bayes, mle, unsure of BK)

	results <- list("mod"=mod, "metab"=data.frame("GPP"=gpp, "R"=resp, "NEP"=nep))
	# attr(results, "lm.mod") <- mod
	return(results)
	
}
