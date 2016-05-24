
# ==================================
# = Temperature Kalman filter/ nll =
# ==================================
KFnllTemp <- function(Params, wtr, watts){
	beta <- Params[1]
	c1 <- Params[2]
	Q <- exp(Params[3])
	H <- exp(Params[4])
	
	alpha <- wtr[1]
	P <- Q
	
	nlls <- rep(0, length(wtr))

	nlls <- kalmanLoopTempR(nlls=nlls, alpha=alpha, c1=c1, P=P, Q=Q, H=H, beta=beta, watts=watts, wtr=wtr)
		
	return(sum(nlls))

}

# ===============================
# = Temperature Kalman Smoother =
# ===============================
KFsmoothTemp <- function(Params, wtr, watts){	
	nobs <- length(wtr)
	t0 <- double(nobs-1)
	
	# Unpack parameters (these were previously fitted)
	beta <- Params[1] # AR1 coefficient (should be less than 1, represents energy loss)
	c1 <- Params[2] # heat gained per watts entering but not leaving layer (accounts for estimate of watts not being perfect)
	Q <- exp(Params[3]) #forced Q to be pos to get rid of warnings; check kfmetab_chla and the NLL function.
	H <- exp(Params[4]) 
	
	# Set first true value equal to first observation
	alpha <- wtr[1]
	
	# Set process covariance, P, equal to Q
	P <- Q # starting value
	
	# Initial values
	aHat <- c(alpha, t0) # aHat[t] == "a[t|t-1]" (estimate of a before updating)
	pHat <- c(P, t0) # pHat[t] == "p[t|t-1]" (estimate of a before updating)
	aVec <- aHat # aVec[t] == "a[t|t]" or "a[t]" (aVec is the "updated" version of aHat)
	pVec <- pHat # pVec[t] == "P[t|t]" or "P[t]" (pVec is the "updated" version of pHat)
	etaVec <- double(nobs)
	
	
	#Pseudocode #3: Set process covariance, P, equal to Q
	P <- Q
	
	for(i in 2:nobs){
		
		#Pseudocode #4a: Predictions
		alpha <- beta*alpha + c1*watts[i-1]
		aHat[i] <- alpha
		
		P <- beta*P*beta + Q
		pHat[i] <- P
		
		#Pseudocode #4b: Update Predictions
		eta <- wtr[i] - alpha
		Eff <- P + H
		alpha <- alpha + P/Eff*eta
		P <- P - P*P/Eff
		
		etaVec[i] <- eta
		aVec[i] <- alpha
		pVec[i] <- P
		}
	
	#Kalman Smoother
	aSmooth <- rep(NA,nobs)
	pSmooth <- rep(NA,nobs)
	aSmooth[nobs] <- aVec[nobs]
	
	for(i in (nobs-1):1){
		pStar <- pVec[i]*beta/pHat[i+1]
		aSmooth[i] <- aVec[i] + pStar*(aSmooth[i+1] - aHat[i+1])
	}
	return(aSmooth) # return smoothed temperature time series
}

# =========================
# = R Function for C loop =
# =========================
kalmanLoopTempR <- function(nlls, alpha, c1,  P, Q, H, beta, watts, wtr){
	nobs <- length(wtr)
	a.loop <- .C("kalmanLoopTempC", nlls=as.double(nlls), as.double(alpha), as.double(c1), as.double(P), as.double(Q), as.double(H), as.double(beta), as.double(watts), as.double(wtr), as.integer(nobs), PACKAGE="LakeMetabolizer")
	return(a.loop[["nlls"]])
}

#'@title Smooth temperature time series using a Kalman filter/ smoother
#'@description Smoothes a temperature time series uses a Kalman filter/ smoother.
#'@usage
#'temp.kalman(wtr, watts, ampH=1, ...)
#'@details 
#'basic model process is \code{x[t] = beta*x[t-1] + c1*watts[t-1]}
#'@param wtr Vector (regular time series) of water temperature in degrees C
#'@param watts estimate of watts entering the layer at each time step, from \link{watts.in}
#'@param ampH factor by which to artificially amplify the observation error variance, H
#'@param ... parameters to be passed to \link{optim}
#'@return a smoothed temperature time series
#'@author Ryan Batt
#'@references 
#'Batt, Ryan D. and Stephen R. Carpenter. 2012. \emph{Free-water lake metabolism: 
#'addressing noisy time series with a Kalman filter}. Limnology and 
#'Oceanography: Methods 10: 20-30. doi: 10.4319/lom.2012.10.20
#'@seealso \link{watts.in} \link{metab.kalman}
#'@export
temp.kalman <- function(wtr, watts, ampH=1, ...){
	# Filter and fit
	guesses <- c(0.9, 1E-4, log(5), log(5))
	fit <- optim(guesses, fn=KFnllTemp, wtr=wtr, watts=watts, ...)
	pars0 <- fit$par
	pars <- c("ar1"=pars0[1], "c1"=pars0[2], "Q"=exp(pars0[3]), "H"=exp(pars0[4])*ampH)
	
	# Smooth
	# smoothTemp <- 
	KFsmoothTemp(pars, wtr=wtr, watts=watts)
		
}


