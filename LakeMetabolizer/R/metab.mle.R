#'@title Metabolism calculated from the maximum likelihood estimates of the parameters in a standard linear regression model
#'@description Process-error-only model with parameters fitted via maximum likelihood estimation (MLE). This function runs the maximum likelihood metabolism model on the supplied gas concentration and other supporting data.
#'@param do.obs Vector of dissolved oxygen concentration observations, \eqn{mg O[2] L^{-1}}{mg O2 / L}
#'@param do.sat Vector of dissolved oxygen saturation values based on water temperature. Calculate using \link{o2.at.sat}
#'@param k.gas Vector of kGAS values calculated from any of the gas flux models 
#'(e.g., \link{k.cole}) and converted to kGAS using \link{k600.2.kGAS}
#'@param z.mix Vector of mixed-layer depths in meters. To calculate, see \link{ts.meta.depths}
#'@param irr Vector of photosynthetically active radiation in \eqn{\mu mol\ m^{-2} s^{-1}}{micro mols / m^2 / s}
#'@param wtr Vector of water temperatures in \eqn{^{\circ}C}{degrees C}. Used in scaling respiration with temperature
#'@param error.type Option specifying if model should assume pure Process Error 'PE' or Observation Error 'OE'. Defaults to observation error 'OE'. 
#'@param ... additional arguments; currently "datetime" is the only recognized argument passed through \code{...}
#'@return
#'A data.frame with columns corresponding to components of metabolism 
#'\describe{
	#'\item{GPP}{numeric estimate of Gross Primary Production, \eqn{mg O_2 L^{-1} d^{-1}}{mg O2 / L / d}}
	#'\item{R}{numeric estimate of Respiration, \eqn{mg O_2 L^{-1} d^{-1}}{mg O2 / L / d}}
	#'\item{NEP}{numeric estimate of Net Ecosystem production, \eqn{mg O_2 L^{-1} d^{-1}}{mg O2 / L / d}}
#'}
#' The maximum likelihood estimates of model parameters can be accessed via \code{attributes(metab.mle(...))[["params"]]}
#' 
#'@details
#'The model has the three parameters, \eqn{c_1, c_2, \epsilon}{c1, c2, epsilon}, and has the form
#'
#'\deqn{v=k.gas/z.mix}{v=k.gas/z.mix}
#'
#'\deqn{a_t = c_1*irr_{t-1} + c_2*log_e(wtr_{t-1}) + v_{t-1}*do.sat_{t-1}}{a[t] = c1*irr[t-1] + c2*log(wtr[t-1]) + v[t-1]*do.sat[t-1]}
#'
#'\deqn{\beta = e^{-v}}{beta = exp(-v)}
#'
#'\deqn{do.obs_t = a_t/v_{t-1} + -e^{-v_{t-1}}*a_t/v_{t-1} + \beta_{t-1}*\do.obs_{t-1} + \epsilon_t}{do.obs[t] = a[t]/v[t-1] + -exp(-v[t-1])*a[t]/v[t-1] + beta[t-1]*do.obs[t-1] + epsilon[t]}
#'
#'
#' The above model is used during model fitting, but if gas flux is not integrated between time steps, those equations simplify to the following:
#'
#' \deqn{F_{t-1} = k.gas_{t-1}*(do.sat_{t-1} - do.obs_{t-1})/z.mix_{t-1}}{F[t-1] = k.gas[t-1]*(do.sat[t-1] - do.obs[t-1])/z.mix[t-1]}
#'
#'\deqn{do.obs_t=do.obs_{t-1}+c_1*irr_{t-1}+c_2*log_e(wtr_{t-1}) + F_{t-1} + \epsilon_t}{do.obs[t] = do.obs[t-1] + c1*irr[t-1] + c2*log(wtr[t-1]) + F[t-1] + epsilon[t]}
#'
#'
#'The parameters are fit using maximum likelihood, and the optimization (minimization of the negative log likelihood function) is performed by \code{optim} using default settings. 
#'
#'GPP is then calculated as \code{mean(c1*irr, na.rm=TRUE)*freq}, where \code{freq} is the number of observations per day, as estimated from the typical size between time steps. Thus, generally \code{freq==length(do.obs)}. 
#'
#'Similarly, R is calculated as \code{mean(c2*log(wtr), na.rm=TRUE)*freq}. 
#'
#'NEP is the sum of GPP and R. 
#'
#'@note Currently, missing values in any arguments will result in an error, so freq must always equal nobs.
#'@author Luke A Winslow, Ryan Batt, GLEON Fellows
#'@references
#'Hanson, PC, SR Carpenter, N Kimura, C Wu, SP Cornelius, TK Kratz. 2008 
#'\emph{Evaluation of metabolism models for free-water dissolved oxygen in lakes}. 
#'Limnology and Oceanography: Methods 6: 454:465.
#'
#'Solomon CT, DA Bruesewitz, DC Richardson, KC Rose, MC Van de Bogert, PC Hanson, TK Kratz, B Larget, 
#'R Adrian, B Leroux Babin, CY Chiu, DP Hamilton, EE Gaiser, S Hendricks, V Istvanovics, A Laas, DM O'Donnell, 
#'ML Pace, E Ryder, PA Staehr, T Torgersen, MJ Vanni, KC Weathers, G Zhuw. 2013. 
#'\emph{Ecosystem Respiration: Drivers of Daily Variability and Background Respiration in Lakes around the Globe}. 
#'Limnology and Oceanography 58 (3): 849:866. doi:10.4319/lo.2013.58.3.0849.
#'
#'@importFrom stats dnorm optim
#'
#'@seealso
#'\link{metab}, \link{metab.bookkeep}, \link{metab.ols}, \link{metab.kalman}, \link{metab.bayesian}
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
#'mod.date = as.POSIXct('2009-07-08', 'GMT')
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
#'metab.mle(doobs[,2], do.sat, k.gas, z.mix[,2], irr[,2], wtr[,3])
#'@export
metab.mle <- function(do.obs, do.sat, k.gas, z.mix, irr, wtr, error.type="OE", ...){

  complete.inputs(do.obs=do.obs, do.sat=do.sat, k.gas=k.gas, 
                  z.mix=z.mix, irr=irr, wtr=wtr, error=TRUE)
  
	match.arg(error.type, choices=c('OE', 'PE'))
	
	nobs <- length(do.obs)
	
	mm.args <- list(...)
  
  if(any(z.mix <= 0)){
    stop("z.mix must be greater than zero.")
  }
	if(any(wtr <= 0)){
		stop("all wtr must be positive.")
	}
  
	if("datetime"%in%names(mm.args)){ # check to see if datetime is in the ... args
		datetime <- mm.args$datetime # extract datetime
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
	
	chk.list <- list(do.obs, irr, do.sat, z.mix, k.gas, wtr)
	if(!all(sapply(chk.list, is.numeric)) || !all(sapply(chk.list, is.vector))){
		stop('All metab.mle inputs must be numeric vectors.')
	}
 
	if(!all(nobs==sapply(chk.list, length))){
		stop('All input data to metab.mle must be the same length')
	}
	
	Q0 <- ((diff(range(do.obs,na.rm=TRUE)) - mean(do.obs,na.rm=TRUE))^2 / length(do.obs))
	
	guesses <- c(1E-4, 1E-4, log(Q0))
	
	#We have a different number of fitted parameters depending on error type of the model
	if(error.type=='OE'){
		guesses <- c(guesses, do.obs[1]) 
	
		fit <- optim(guesses, fn=mleNllOE, do.obs=do.obs, do.sat=do.sat, k.gas=(k.gas/freq), z.mix=z.mix, irr=irr, wtr=wtr)
		
		pars0 <- fit$par
		pars <- c("gppCoeff"=pars0[1], "rCoeff"=pars0[2], "Q"=exp(pars0[3]), "nll"=fit$value, "doInit"=pars0[4])
		
	}else if(error.type=='PE'){
		guesses <- c(guesses) 
		
		fit <- optim(guesses, fn=mleNllPE, do.obs=do.obs, do.sat=do.sat, k.gas=(k.gas/freq), z.mix=z.mix, irr=irr, wtr=wtr)
		
		pars0 <- fit$par
		pars <- c("gppCoeff"=pars0[1], "rCoeff"=pars0[2], "Q"=exp(pars0[3]), "nll"=fit$value)
		
	}else{
		stop("error.type must be either 'OE' or 'PE', Observation Error or Process Error respectively.")
	}
	
	# ====================================
	# = Use fits to calculate metabolism =
	# ====================================
	GPP <- mean(pars[1]*irr, na.rm=TRUE) * freq
	R <- mean(pars[2]*log(wtr), na.rm=TRUE) * freq
	
	return(list("params"=pars, "metab"=c("GPP"=GPP,"R"=R,"NEP"=GPP+R)))
}

# ============================================
# = The R loop for Observation Error mle NLL =
# ============================================
mleLoopOE <- function(alpha, doobs, c1, c2, beta, irr, wtr, kz, dosat){
	nobs <- length(doobs)
	a.loop <- .C("mleLoopCoe", alpha=as.double(alpha), as.double(doobs), as.double(c1), as.double(c2), as.double(beta), as.double(irr), as.double(wtr), as.double(kz), as.double(dosat), as.integer(nobs), PACKAGE="LakeMetabolizer")
	return(a.loop[["alpha"]])
}

# ============================================
# = The R loop for Process Error mle NLL =
# ============================================
mleLoopPE <- function(alpha, doobs, c1, c2, beta, irr, wtr, kz, dosat){
	nobs <- length(doobs)
	a.loop <- .C("mleLoopCpe", alpha=as.double(alpha), as.double(doobs), as.double(c1), as.double(c2), as.double(beta), as.double(irr), as.double(wtr), as.double(kz), as.double(dosat), as.integer(nobs), PACKAGE="LakeMetabolizer")
	return(a.loop[["alpha"]])
}

# ====================
# = mle NLL function =
# ====================
mleNllPE <- function(Params, do.obs, do.sat, k.gas, z.mix, irr, wtr){
	c1 <- Params[1] #PAR coeff
	c2 <- Params[2] #log(Temp) coeff
	Q <- exp(Params[3]) # Variance of the process error
	
	# See KalmanDO_smooth.R comments for explanation of beta
	kz <- k.gas/z.mix # K and Zmix are both vector of length nobs
	beta <- exp(-kz) # This beta is for using the differential equation form
	
	# Set first true value equal to first observation
	alpha <- rep(0, length(do.obs))
	alpha[1] <- do.obs[1]#Let's give this model some starting values

	#R version of C loop
	#for(i in 2:length(do.obs)){
	#	a1 <- c1*irr[i-1] + c2*log(wtr[i-1]) + kz[i-1]*do.sat[i-1]
	#	alpha[i] <- a1/kz[i-1] + -exp(-kz[i-1])*a1/kz[i-1] + beta[i-1]*alpha[i-1] # NOTE: beta==exp(-kz); kz=K/Zmix
	#}
	alpha <- mleLoopPE(alpha=alpha, doobs=do.obs, c1=c1, c2=c2, beta=beta, irr=irr, wtr=wtr, kz=kz, dosat=do.sat)
	
	return(-sum(dnorm(do.obs, alpha, sd=sqrt(Q), log=TRUE), na.rm=TRUE))
}#End function

# ====================
# = mle NLL function =
# ====================
mleNllOE <- function(Params, do.obs, do.sat, k.gas, z.mix, irr, wtr, error.type){
	c1 <- Params[1] #PAR coeff
	c2 <- Params[2] #log(Temp) coeff
	Q <- exp(Params[3]) # Variance of the process error
	
	# See KalmanDO_smooth.R comments for explanation of beta
	kz <- k.gas/z.mix # K and Zmix are both vector of length nobs
	beta <- exp(-kz) # This beta is for using the differential equation form
	
	# Set first true value equal to first observation
	alpha <- rep(0, length(do.obs))
	alpha[1] <- Params[4] #Free varying initial DO value
	
	alpha <- mleLoopOE(alpha=alpha, doobs=do.obs, c1=c1, c2=c2, beta=beta, irr=irr, wtr=wtr, kz=kz, dosat=do.sat)
	
	return(-sum(dnorm(do.obs, alpha, sd=sqrt(Q), log=TRUE), na.rm=TRUE))
}#End function
