#'@title Metabolism calculated from parameters estimated using a Kalman filter
#'@description A state space model accounting for process and observation error, with the maximum likelihood of parameteres estimated using a Kalman filter.
#'Also provides a smoothed time series of oxygen concentration.
#'@param do.obs Vector of dissovled oxygen concentration observations, \eqn{mg O[2] L^{-1}}{mg O2 / L}
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
#'
#' Use \link{attributes} to access more model output:
#'\item{smoothDO}{smoothed time series of oxygen concentration (\eqn{mg O[2] L^{-1}}{mg O2 / L}), from Kalman smoother}
#'\item{params}{parameters estimated by the Kalman filter (\eqn{c_1, c_2, Q, H}{c1, c2, Q, H})}
#' 
#'@details
#'The model has four parameters, \eqn{c_1, c_2, Q, H}{c1, c2, Q, H}, and consists of equations involving the prediction of upcoming state conditional on information of the previous state (\eqn{a_{t|t-1}}{a[t|t-1]}, \eqn{P_{t|t-1}}{P[t|t-1]}), as well as updates of those predictions that are conditional upon information of the current state (\eqn{a_{t|t}}{a[t|t]}, \eqn{P_{t|t}}{P[t|t]}). \eqn{a} is the 
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
#'@references
#'Batt, Ryan D. and Stephen R. Carpenter. 2012. \emph{Free-water lake metabolism: 
#'addressing noisy time series with a Kalman filter}. Limnology and Oceanography: Methods 10: 20-30. doi: 10.4319/lom.2012.10.20
#'@seealso
#'\link{temp.kalman}, \link{watts.in}, \link{metab}, \link{metab.bookkeep}, \link{metab.ols}, \link{metab.mle}, \link{metab.bayesian}
#'@author Ryan Batt, Luke A. Winslow
#'@note If observation error is substantial, consider applying a Kalman filter to the water temperature time series by supplying 
#' \code{wtr} as the output from \link{temp.kalman}
#'@examples
#'library(rLakeAnalyzer)
#'doobs <- load.ts(system.file('extdata', 
#'                             'sparkling.doobs', package="LakeMetabolizer"))
#'wtr <- load.ts(system.file('extdata', 
#'                           'sparkling.wtr', package="LakeMetabolizer"))
#'wnd <- load.ts(system.file('extdata', 
#'                           'sparkling.wnd', package="LakeMetabolizer"))
#'irr <- load.ts(system.file('extdata', 
#'                           'sparkling.par', package="LakeMetabolizer"))
#'
#'#Subset a day
#'Sys.setenv(TZ='GMT')
#'mod.date <- as.POSIXct('2009-07-08', 'GMT')
#'doobs <- doobs[trunc(doobs$datetime, 'day') == mod.date, ]
#'wtr <- wtr[trunc(wtr$datetime, 'day') == mod.date, ]
#'wnd <- wnd[trunc(wnd$datetime, 'day') == mod.date, ]
#'irr <- irr[trunc(irr$datetime, 'day') == mod.date, ]
#'
#'k600 <- k.cole.base(wnd[,2])
#'k.gas <- k600.2.kGAS.base(k600, wtr[,3], 'O2')
#'do.sat <- o2.at.sat.base(wtr[,3], altitude=300)
#'
#'metab.kalman(irr=irr[,2], z.mix=rep(1, length(k.gas)), 
#'             do.sat=do.sat, wtr=wtr[,2],
#'             k.gas=k.gas, do.obs=doobs[,2])
#'@export
metab.kalman <- function(do.obs, do.sat, k.gas, z.mix, irr, wtr, ...){
  
  complete.inputs(do.obs=do.obs, do.sat=do.sat, k.gas=k.gas, 
                  z.mix=z.mix, irr=irr, wtr=wtr, error=TRUE)
  
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
  
  # Filter and fit
  guesses <- c(1E-4,1E-4,log(5),log(5))
  fit <- optim(guesses, fn=KFnllDO, do.obs=do.obs, do.sat=do.sat, k.gas=(k.gas/freq), z.mix=z.mix, irr=irr, wtr=wtr)
  pars0 <- fit$par
  pars <- c("gppCoeff"=pars0[1], "rCoeff"=pars0[2], "Q"=exp(pars0[3]), "H"=exp(pars0[4]))
  
  # Smooth
  smoothDO <- KFsmoothDO(pars, do.obs=do.obs, do.sat=do.sat, k.gas=(k.gas/freq), z.mix=z.mix, irr=irr, wtr=wtr)
  
  # Use fits to calculate metabolism
  
  GPP <- mean(pars[1]*irr, na.rm=TRUE) * freq
  R <- mean(pars[2]*log(wtr), na.rm=TRUE) * freq
  
  return(list("smoothDO"=smoothDO,"params"=pars, "metab"=c("GPP"=GPP,"R"=R, "NEP"=GPP+R)))
}




# ======================
# = Kalman filter/ nll =
# ======================
# Main recursion written in C
KFnllDO <- function(Params, do.obs, do.sat, k.gas, z.mix, irr, wtr){
	
	# ===========================
	# = Unpack and set initials =
	# ===========================
	#!Pseudocode #1: Initial guesses for B, C, and Q t
	c1 <- Params[1] #PAR coeff
	c2 <- Params[2] #log(Temp) coeff
	Q <- exp(Params[3]) # Variance of the process error
	H <- exp(Params[4]) # Variance of observation error
	
	# See KalmanDO_smooth.R comments for explanation of beta
	kz <- k.gas/z.mix # K and Zmix are both vector of length nobs
	# beta <- 1-kz # beta is a vector of length nobs (this beta is for difference equation form)
	beta <- exp(-kz) # This beta is for using the differential equation form
	
	# Set first true value equal to first observation
	alpha <- do.obs[1]#Let's give this model some starting values
	
	# Set process covariance, P, equal to Q
	P <- Q #starting value
	
	# Empty vector for nll's
	nlls <- rep(0,length(do.obs))
		
	# ==================
	# = Main Recursion =
	# ==================
	# ===============
	# = Predictions =
	# ===============
	# Equations for Predictions from Harvey
	# a[t|t-1] = T[t]*a[t-1] + c[t] Harvey pg 105 eq. 3.2.2a
	# P[t|t-1] = T[t]*P[t-1]*T'[t] + R[t]*Q[t]*R'[t] Harvey pg 106 eq. 3.2.2b
	
	# Predictions where gas flux not split into beta etc. (I'm pretty sure this is wrong):
	# Uk <- K[i-1]*(do.sat[i-1] - alpha)/Zmix[i-1]
	# alpha <- alpha + c1*irr[i-1] + c2*log(wtr[i-1]) + Uk
	# P <- (Uk*P*Uk) + Q
	
	# Predictions where gas flux is split into beta:
	# Difference Equation Version (original):
	# alpha <- beta[i-1]*alpha + c1*irr[i-1] + c2*log(wtr[i-1]) + kz[i-1]*do.sat[i-1]
	
	# Differential Equation Version:
	# Gordon's code (from BayesMetabSS_indivProcessErr_072012.txt):
	# alpha[f,j] <- P[f,j] - rho[f,j] + KO2[f,j]* DOSat[f,j];
	# DOHat[f+1,j] <- (alpha[f,j]-(alpha[f,j]-KO2[f,j]*DOTrue[f,j])*exp(-KO2[f,j]))/KO2[f,j];
	
	# Define Gordon's "alpha" as "a1":
	# a1 = c1*irr + c2*log(wtr) + kz*do.sat # my a1 is his alpha
	
	# Coefficients in front of alpha are defined above as "beta", and play a particular role in propagating uncertainty
	# Need to algebraically regarrange Gordon's equation for DOHat (my "alpha" is his "DOHat")
	# Note that alpha and a1 are rewritten each iteration of the loop...
	# Therefore, alpha on the left side is "alpha[t]", and alpha on the right side is "alpha[t-1]"
	
	# alpha = (a1 - (a1 - kz*alpha)*exp(-kz))/kz # my alpha is his DOHat
	# alpha = (a1 + (-a1 + kz*alpha)*exp(-kz))/kz # redistribute -1
	# alpha = (a1 + exp(-kz)*(-a1 + kz*alpha))/kz # regarrange
	# alpha = (a1 + -exp(-kz)*a1 + exp(-kz)*kz*alpha)/kz # multiply "exp(kz)" through
	# alpha = a1/kz + -exp(-kz)*a1/kz + exp(-kz)*alpha # multiply "/kz" through
	
	# ======================
	# = Updating Equations =
	# ======================
	# Updating Equations from Harvey
	# a[t] = a[t|t-1] + P[t|t-1]*Z'[t]*F[t]^-1(y[t] - Z[t]*a[t|t-1] - d[t]). Harvey, page 106, 3.2.3a
	# P[t] = P[t|t-1] - P[t|t-1]*Z'[t]*F[t]^-1*Z[t]*P[t|t-1] Harvey, page 106, eq. 3.2.3b
	# F[t] = Z[t]*P[t|t-1]*Z'[t] + H[t] Harvey, page 106, eq. 3.2.3c
	# kalmanLoopC(double *alpha, double *doobs, double *c1, double *c2, double *P, double *Q, double *H,  double *beta, double *irr, double *wtr, double *kz, double *dosat, int *nobs)
	nlls <- kalmanLoopR(nlls=nlls, alpha=alpha, doobs=do.obs, c1=c1, c2=c2, P=P, Q=Q, H=H, beta=beta, irr=irr, wtr=wtr, kz=kz, dosat=do.sat)

	return(sum(nlls)) # return the sum of nll's
}#End function



# ===================
# = Kalman Smoother =
# ===================
KFsmoothDO <- function(Params, do.obs, do.sat, k.gas, z.mix, irr, wtr, Hfac=NULL){
	nobs <- length(do.obs)
	d0 <- double(nobs-1)
	# beta <- 1-KO2zmix #do.obs_t = 1*do.obs_t-1 + -KO2zmix*do.obs_t-1 + Sea%*%Ewe + eta === (1-KO2zmix)*do.obs_t-1.....
	
	# Unpack parameters (these were previously fitted)
	c1 <- Params[1] # irr Coeff
	c2 <- Params[2] # log(wtr) Coeff
	Q <- Params[3] # Variance of the Process Error
	if(is.null(Hfac)){
		H <- Params[4]
	}else{
		H <- Params[4]*Hfac
	}
	 # Variance of Observation Error
	
	# Need to define portion of K multiplied by state variable (DO)
	# Gas flux = K[t-1](do.sat[t-1] - alpha[t-1])/Zmix[t-1]
	# Gas flux = (K[t-1]/Zmix[t-1])*(do.sat[t-1]) - (K[t-1]/Zmix[t-1])*(alpha[t-1])
	# Note that K/Z is essentially a coefficient sitting in front of alpha, the estimate of DO
	# Therefore,
	# alpha[t] = 1*alpha[t-1] + c1*irr[t-1] + c2*log(wtr[t-1]) + K[t-1](do.sat[t-1] - alpha[t-1])/Zmix[t-1]
	# Becomes:
	# kz[t] = k[t]/Zmix[t]
	# alpha[t] = 1*alpha[t-1] + -kz[t-1]*alpha[t-1] + c1*irr[t-1] + c2*log(wtr[t-1]) + kz[t-1]*do.sat[t-1]
	# Or,
	# alpha[t] = (1-kz[t-1])*alpha[t-1] + c1*irr[t-1] + c2*log(wtr[t-1]) + kz[t-1]*do.sat[t-1]
	# Defining kz and redefining (1-kz[t]) as beta[t]:
	kz <- k.gas/z.mix # K and Zmix are both vector of length nobs
	# beta <- 1-kz # beta is a vector of length nobs (this beta is for difference equation form)
	beta <- exp(-kz) # This beta is for using the differential equation form
	
	# Set first true value equal to first observation
	alpha <- do.obs[1]
	
	# Set process covariance, P, equal to Q
	P <- Q # starting value
	
	# Initial values
	aHat <- c(alpha, d0) # aHat[t] == "a[t|t-1]" (estimate of a before updating)
	pHat <- c(P, d0) # pHat[t] == "p[t|t-1]" (estimate of a before updating)
	aVec <- aHat # aVec[t] == "a[t|t]" or "a[t]" (aVec is the "updated" version of aHat)
	pVec <- pHat # pVec[t] == "P[t|t]" or "P[t]" (pVec is the "updated" version of pHat)
	etaVec <- double(nobs)
	
	for(i in 2:nobs){
		# ===============
		# = Predictions =
		# ===============
		# Equations for Predictions from Harvey
		# a[t|t-1] = T[t]*a[t-1] + c[t] Harvey pg 105 eq. 3.2.2a
		# P[t|t-1] = T[t]*P[t-1]*T'[t] + R[t]*Q[t]*R'[t] Harvey pg 106 eq. 3.2.2b
		
		# Predictions where gas flux not split into beta etc.:
		# Uk <- K[i-1]*(do.sat[i-1] - alpha)/Zmix[i-1]
		# alpha <- alpha + c1*irr[i-1] + c2*log(wtr[i-1]) + Uk
		# aHat[i] <- alpha
		# P <- (Uk*P*Uk) + Q
		# pHat[i] <- P
		
		# Predictions where gas flux is split into beta (see explanation above):
		
		# Difference Equation Version:
		# alpha <- beta[i-1]*alpha + c1*irr[i-1] + c2*log(wtr[i-1]) + kz[i-1]*do.sat[i-1]
		
		# Differential Equation Version (see kalmanDO_nll.R for explanation):
		if(is.finite(1/kz[i-1])){
			
			a1 <- c1*irr[i-1] + c2*log(wtr[i-1]) + kz[i-1]*do.sat[i-1]	
			alpha <- a1/kz[i-1] + -beta[i-1]*a1/kz[i-1] + beta[i-1]*alpha # NOTE: beta==exp(-kz); kz=K/Zmix
			
		}else{
			
			alpha <- c1*irr[i-1] + c2*log(wtr[i-1])
			
		}

		
		aHat[i] <- alpha
		P <- (beta[i-1]*P*beta[i-1]) + Q
		pHat[i] <- P
	
		# ======================
		# = Updating Equations =
		# ======================
		# Updating Equations from Harvey
		# a[t] = a[t|t-1] + P[t|t-1]*Z'[t]*F[t]^-1(y[t] - Z[t]*a[t|t-1] - d[t]). Harvey, page 106, 3.2.3a
		# P[t] = P[t|t-1] - P[t|t-1]*Z'[t]*F[t]^-1*Z[t]*P[t|t-1] Harvey, page 106, eq. 3.2.3b
		# F[t] = Z[t]*P[t|t-1]*Z'[t] + H[t] Harvey, page 106, eq. 3.2.3c
				
		eta <- do.obs[i] - alpha
		Eff <- P + H
		alpha <- alpha + P/Eff*eta
		P <- P - P*P/Eff

		aVec[i] <- alpha
		pVec[i] <- P
		etaVec[i] <- eta
	}
	
	#Kalman Smoother
	aSmooth <- rep(NA,nobs)
	Psmooth <- rep(NA,nobs)
	aSmooth[nobs] <- aVec[nobs] # "starting" value for smoother (smoother starts at end and works backwards)
	# pSmooth[nobs] <- pVec[nobs]
	
	# Filtering is informed by past information
	# Smoothing includes the information from filtering (estimates of parameters), but also future information.
	# "The aim of filtering is to find the expected value of the state vector, alpha[t], conditional on the information available at time t, that is E(alpha[t]|Y[t]). The aim of smoothing is to take account of the information made available after time t. The mean of the distribution of alpha[t], conditional on all the sample, may be written as E(alpha[t]|Y[T]) and is known as the smoothed estimate. THe corresponding estimator is called the SMOOTHER. Since the smoother is based on more information than the filtered estimator, it will have a MSE which, in general, is smaller than that of the filtered estimator; it cannot be greater." ~ Harvey 1989, pgs 149-150.
	 #a[t|T] = a[t] + Pstar[t]*(a[t+1|T] - T[t+1]*a[t])
	# P[t|T] = P[t] + Pstar[t]*(P[t+1|T] - P[t+1|t])*Pstar[t]
	# Pstar[t] = P[t]*T[t+1]/P[t+1|t]
	# t is current time step, T is last time step (when in []), T contains AR parameters (when NOT in [])
	
	for(i in length(d0):1){
		pStar <- pVec[i]*beta[i+1]/pHat[i+1]
		aSmooth[i] <- aVec[i] + pStar*(aSmooth[i+1] - aHat[i+1])
		
		# CAN ALSO SMOOTH P, WHICH GIVES THE SMOOTHED COVARIANCE MATRIX (not a matrix for univariate; gives estimate of accuracy of state estimate)
		# pSmooth[i] <- pVec[i] + pStar*(pSmooth[i+1] - pHat[i+1])*pStar
		}
	return(aSmooth)	# return smoothed DO time series
}

# ===========================================
# = R function that calls the C loop for KF =
# ===========================================
# kalmanLoopC(double *alpha, double *doobs, double *c1, double *c2, double *P, double *Q, double *H,  double *beta, double *irr, double *wtr, double *kz, double *dosat, int *nobs)
kalmanLoopR <- function(nlls, alpha, doobs, c1, c2, P, Q, H, beta, irr, wtr, kz, dosat){
	nobs <- length(doobs)
	a.loop <- .C("kalmanLoopC", nlls=as.double(nlls), as.double(alpha), as.double(doobs), as.double(c1), as.double(c2), as.double(P), as.double(Q), as.double(H), as.double(beta), as.double(irr), as.double(wtr), as.double(kz), as.double(dosat), as.integer(nobs), PACKAGE="LakeMetabolizer")
	return(a.loop[["nlls"]])
}

