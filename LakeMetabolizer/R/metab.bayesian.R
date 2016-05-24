
if(getRversion() >= "2.15.1")  utils::globalVariables("jags")

# ====================================
# = Function to write the jags model =
# ====================================
bayes.makeModel <- function(k.gas){
  
	if(!requireNamespace("R2jags")){
    stop("metab.bayesian requires R2jags.\ninstall.packages('R2jags')\n #Also install JAGS (http://mcmc-jags.sourceforge.net/)")
	}
	
	finite.1oK <- is.finite(1/k.gas)
	
	choose.allK <- all(finite.1oK)
	choose.noK <- all(!finite.1oK)
	choose.bothK <- any(finite.1oK) & any(!finite.1oK)
	
	choice.mod <- c("allK", "noK", "bothK")[c(choose.allK, choose.noK, choose.bothK)]

	# Write the appropriate bayesian model into a temporary file
	#modfile <- tempfile('jags.metab.bayes')
	jags.dir = system.file('jags', package='LakeMetabolizer')
	switch(choice.mod,
		allK = {modfile = file.path(jags.dir, 'bayes.mod.allk.JAGS'); print("using allK model")},
		noK = {modfile = file.path(jags.dir, 'bayes.mod.nok.JAGS'); print("using noK model")},
		bothK = {modfile = file.path(jags.dir, 'bayes.mod.bothk.JAGS'); print("using bothK model")}
	)

	return(modfile)
}
# 
# 
# # ==================================
# # = Bayes model for all non-zero K =
# # ==================================
# bayes.mod.allK <- function(){
# 	
# 	# model process
# 	for(i in 2:N){	
# 		Y[i] ~ dnorm(a[i], tauV) # observations (Y) are distributed with mean equivalent to true values, and precision tauV (1/tauV is variance of observation error)
# 		K[i-1] ~ dnorm(kP[i-1, 1], 1/kP[i-1, 2]) #distributin on K
# 		
# 		kz[i-1] <- K[i-1]/Zmix[i-1]
# 		
# 		a1[i] <- U[i-1,]%*%C + kz[i-1]*satO[i-1]
# 		aHat[i] <- a1[i]/kz[i-1] + -exp(-1*kz[i-1])*a1[i]/kz[i-1] + exp(-kz[i-1])*a[i-1]
# 		
# 		a[i] ~ dnorm(aHat[i], tauW) # true values have a mean equivalent to estimated values, but accompanied by process error (process precision is tauW)
# 	}
# 	
# 	# Starting values
# 	a[1] <- a0
# 	
# 	#Priors on regression coefficients
# 	C[1] ~ dnorm(cP[1,1], 1/cP[1,2])
# 	C[2] ~ dnorm(cP[2,1], 1/cP[2,2])
# 	
# 	#Prior on errors
# 	tauV ~ dgamma(1.0E-3, 1.0E-3)
# 	tauW ~ dgamma(1.0E-3, 1.0E-3)
# 	sigmaV <- 1/sqrt(tauV)
# 	sigmaW <- 1/sqrt(tauW)	
# }
# 
# # ===============================
# # = Bayes model w/ K turned off =
# # ===============================
# bayes.mod.noK <- function(){
# 
# 	# model process
# 	for(i in 2:N){	
# 		Y[i] ~ dnorm(a[i], tauV) # observations (Y) are distributed with mean equivalent to true values, and precision tauV (1/tauV is variance of observation error)
# 		# NOT USED
# 		K[i-1] ~ dnorm(kP[i-1, 1], 1/kP[i-1, 2]) #distribution on K
# 		
# 		aHat[i] <- a[i-1] + U[i-1,]%*%C # the process
# 
# 		a[i] ~ dnorm(aHat[i], tauW) # true values have a mean equivalent to estimated values, but accompanied by process error (process precision is tauW)
# 	}
# 
# 	# Starting values
# 	a[1] <- a0
# 
# 	#Priors on regression coefficients
# 	C[1] ~ dnorm(cP[1,1], 1/cP[1,2])
# 	C[2] ~ dnorm(cP[2,1], 1/cP[2,2])
# 
# 
# 	#Prior on errors
# 	tauV ~ dgamma(1.0E-3, 1.0E-3)
# 	tauW ~ dgamma(1.0E-3, 1.0E-3)
# 	sigmaV <- 1/sqrt(tauV)
# 	sigmaW <- 1/sqrt(tauW)	
# }
# 
# # ===============================================
# # = Bayes model that can handle both 0 and !0 K =
# # ===============================================
# bayes.mod.bothK <- function(){
# 	
# 	# model process
# 	for(i in 2:N){	
# 		Y[i] ~ dnorm(a[i], tauV) # observations (Y) are distributed with mean equivalent to true values, and precision tauV (1/tauV is variance of observation error)
# 		K[i-1] ~ dnorm(kP[i-1, 1], 1/kP[i-1, 2]) #distributin on K
# 		
# 		# jags cannot handle an expression that includes division by 0 (so can't do blah <- ifelse(kz==0, 1, 1/kz), b/c it'll do 1/kz even when kz==0)
# 		# so to avoid division by 0 when k.gas==0, have to make a "safe" kz that is 1 if kz is 0
# 		# when kzSafe is 1 (to avoid division by 0), that means that we just want aHat to be the bio process
# 		# so we have to cancel out all math that is done by a kzSafe==1 (b/c it is bogus) by multiplying by 0, and add bio process
# 		# but if kzSafe!=0, we should multiply the math involving kzSafe by 1 (to keep it), and multiply the added bio process by 0 (to remove it)
# 		# this is a pretty hacked solution, but I don't see another choice (the only control flow in jags is ifelse, no if(){}else{})
# 			
# 		kz[i-1] <- K[i-1]/Zmix[i-1]
# 		kzSafe[i-1] <- ifelse(kz[i-1]==0, 1, kz[i-1]) # if kz is 0, change to 1
# 		kzCancel[i-1] <- ifelse(kzSafe[i-1]==1, 0, 1) # if kzSafe is 1 (meaning kz is 0), kzCancel needs to be 0
# 		pCancel[i-1] <- ifelse(kzSafe[i-1]==1, 1, 0) # if kzSafe is just kz, then math involving kzSafe is not bogus, and need to cancel the added process
# 		
# 		a1[i] <- U[i-1,]%*%C + kz[i-1]*satO[i-1]
# 		aHat[i] <- (a1[i]/kzSafe[i-1] + -exp(-1*kzSafe[i-1])*a1[i]/kzSafe[i-1] + exp(-kzSafe[i-1])*a[i-1])*kzCancel[i-1] + (a[i-1] + a1[i])*pCancel[i-1]
# 
# 		a[i] ~ dnorm(aHat[i], tauW) # true values have a mean equivalent to estimated values, but accompanied by process error (process precision is tauW)
# 	}
# 	
# 
# 	# Starting values
# 	a[1] <- a0
# 	
# 	#Priors on regression coefficients
# 	C[1] ~ dnorm(cP[1,1], 1/cP[1,2])
# 	C[2] ~ dnorm(cP[2,1], 1/cP[2,2])
# 
# 	#Prior on errors
# 	tauV ~ dgamma(1.0E-3, 1.0E-3)
# 	tauW ~ dgamma(1.0E-3, 1.0E-3)
# 	sigmaV <- 1/sqrt(tauV)
# 	sigmaW <- 1/sqrt(tauW)	
# }


# ================================
# = Supply Data and run bayesFit =
# ================================
bayesFit <- function(data, params, mf, tend="median", ...){ #function that writes jags model, traces params, supplies data, etc
	
	bf.args <- list(...)
	
	jags.m <- R2jags::jags(data, NULL, parameters.to.save=params, mf)

	tF <- function(x, tend){ # tendency function
		switch(tend,
			median=median(x), #median
			mean=mean(x), #mean
			mode1 = unique(x)[which.max(tabulate(match(x, unique(x))))], # mode --- most frequently observed value
			mode2 = { # mode --- based on the highest peak of the posterior probability (from density plot)
				xd <- density(x)
				xd$x[which.max(xd$y)]
			}
			)
	}
	# medSim <- matrix(apply(jags2.m$BUGSoutput$sims.matrix, 2, median)[-(1)], nrow=115, ncol=3)
	simOut <- jags.m$BUGSoutput$sims.matrix
	ctSim <- apply(simOut, 2, tF, tend) # the central tendency metric
	sdSim <- apply(simOut, 2, sd)

	#Figure out the order of the sims.matrix columns ...
	n.obs <- length(data$U[,1])
	GPP <- mean(ctSim[1]*data$U[,1], na.rm=TRUE) * n.obs # gpp coef * par, then sum
	R <- mean(ctSim[2]*data$U[,2], na.rm=TRUE) * n.obs # r coef * log(temp), then sum

	GPPsd <- sqrt(sum(sdSim[1]^2*data$U[,1]^2))
	Rsd <- sqrt(sum(sdSim[2]^2*data$U[,2]^2))
	NEPsd <- sqrt(GPPsd^2 + Rsd^2)

	return(list(
		"model" = jags.m, 
		"params" = ctSim[1:2], 
		"metab.sd" = matrix(c(GPPsd, Rsd, NEPsd), nrow=1, dimnames=list(NULL, c("GPPsd", "Rsd", "NEPsd"))),
		"metab" = matrix(c(GPP, R, GPP+R), nrow=1, dimnames=list(NULL, c("GPP", "R", "NEP")))
	)) # need to clean up format, and maybe include a return of the sd's of the estimates
}

#'@title Metabolism model based on a bayesian parameter estimation framework
#'@description This function runs the bayesian metabolism model on the supplied 
#'gas concentration and other supporting data. This allows for both estimates of 
#'metabolism along with uncertainty around the parameters.
#'@param do.obs Vector of dissovled oxygen concentration observations, mg L^-1
#'@param do.sat Vector of dissolved oxygen saturation values based on water temperature. Calculate using \link{o2.at.sat}
#'@param k.gas Vector of kGAS values calculated from any of the gas flux models 
#'(e.g., \link{k.cole}) and converted to kGAS using \link{k600.2.kGAS}
#'@param z.mix Vector of mixed-layer depths in meters. To calculate, see \link{ts.meta.depths}
#'@param irr Vector of photosynthetically active radiation in \eqn{\mu mol\ m^{-2} s^{-1}}{micro mols / m^2 / s}
#'@param wtr Vector of water temperatures in \eqn{^{\circ}C}{degrees C}. Used in scaling respiration with temperature
#'@param priors Parameter priors supplied as a named numeric vector (example: c("gppMu"=0, "gppSig2"=1E5, "rMu"=0, "rSig2"=1E5, "kSig2"=NA))
#'@param ... additional arguments; currently "datetime" is the only recognized argument passed through \code{...}
#' @return
#' A list of length 4 with components:
	#' \item{model}{the jags model, including posterior draws (see \link[R2jags]{jags})}
	#' \item{params}{parameter estimates of interest from model (medians)}
	#' \item{metab.sd}{standard deviation of metabolism estimates}
	#' \item{metab}{daily metabolism estimates as a data.frame with columns corresponding to
	#' \describe{
		#' \item{\code{GPP}}{numeric estimate of Gross Primary Production, \eqn{mg O_2 L^{-1} d^{-1}}{mg O2 / L / d}}  
		#' \item{\code{R}}{numeric estimate of Respiration, \eqn{mg O_2 L^{-1} d^{-1}}{mg O2 / L / d}}  
		#' \item{\code{NEP}}{numeric estimate of Net Ecosystem production, \eqn{mg O_2 L^{-1} d^{-1}}{mg O2 / L / d}}  
	#' }}


#'@references
#'Holtgrieve, Gordon W., Daniel E. Schindler, Trevor a. Branch, and Z. Teresa A'mar. 
#'2010. \emph{Simultaneous Quantification of Aquatic Ecosystem Metabolism and Reaeration 
#'Using a Bayesian Statistical Model of Oxygen Dynamics}. 
#'Limnology and Oceanography 55 (3): 1047-1062. doi:10.4319/lo.2010.55.3.1047. 
#'http://www.aslo.org/lo/toc/vol_55/issue_3/1047.html.
#'@seealso
#'\link{metab.mle}, \link{metab.bookkeep}, \link{metab.kalman}
#'@author Ryan Batt, Luke A. Winslow
#'
#'@importFrom stats sd density median
#'
#'@examples
#'\dontrun{
#'library(rLakeAnalyzer)
#'
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
#'
#'k600 = k.cole.base(wnd[,2])
#'k.gas = k600.2.kGAS.base(k600, wtr[,3], 'O2')
#'do.sat = o2.at.sat(wtr[,1:2], altitude=300)
#'
#'metab.bayesian(irr=irr[,2], z.mix=rep(1, length(k.gas)), 
#'               do.sat=do.sat[,2], wtr=wtr[,2],
#'               k.gas=k.gas, do.obs=doobs[,2])
#'}
#'@export
metab.bayesian <- function(do.obs, do.sat, k.gas, z.mix, irr, wtr, priors, ...){
	
  complete.inputs(do.obs=do.obs, do.sat=do.sat, k.gas=k.gas, 
                  z.mix=z.mix, irr=irr, wtr=wtr, error=TRUE)
  
  if(any(z.mix <= 0)){
		stop("z.mix must be greater than zero.")
	}
	if(any(wtr <= 0)){
		stop("all wtr must be positive.")
	}
	
	mb.args <- list(...)
	nobs <- length(do.obs)
	# =========================================
	# = Check for datetime and claculate freq =
	# =========================================
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
	
	# ======================================
	# = # Check for priors supplied in ... =
	# ======================================
	if(!missing(priors)){
		t.priors <- priors
		p.name.logic <- all(c(c("gppMu", "gppSig2", "rMu", "rSig2", "kSig2"))%in%names(t.priors))
		p.class.logic <- is.integer(t.priors) | is.numeric(t.priors)
		if(p.name.logic & p.class.logic){
			priors <- t.priors
		}else{
			stop("supplied priors was not a (properly) named numeric/integer vector")
		}
	}else{
		priors <- c("gppMu"=0, "gppSig2"=1E5, "rMu"=0, "rSig2"=1E5, "kSig2"=NA)
	}
		
	if(!all(c(is.numeric(do.obs), is.numeric(do.sat), is.numeric(k.gas), is.numeric(z.mix), is.numeric(irr), is.numeric(wtr)))){
		stop('All inputs to metab.bayes must be numeric vectors')
	}
	
  requireNamespace("R2jags")
	
	# Define model and write to file
	# Model choice depends on k values (all 0, all non-0, mixture)
	modfile <- bayes.makeModel(k.gas=(k.gas/freq))
	
	# ===========================================
	# = Define objects to be used in jags model =
	# ===========================================
	#Supply elements of U (PAR, log(temp))
	U <- matrix(NA, nrow=length(irr), ncol=2)
	U[,1] <- irr # PAR Values
	U[,2] <- log(wtr) # log(temp) values

	# Priors (kP) for k.gas
	kP <- matrix(NA, nrow=length(k.gas), ncol=2)
	# variances for K
	kP[,1] <- k.gas/freq # means for K, adjusted to sampling frequency (velocity denominator in time step, not day)
	if(is.na(priors["kSig2"])){
		k0.logic <- !is.finite(1/kP[,1]) # test for when k is 0
		kP[,2] <- sum(kP[,1])/sum(!k0.logic)*0.1 # k variance = mean of the non-zero K, times 0.1
		kP[k0.logic,2] <- 1E-9
	}else{
		kP[,2] <- priors["kSig2"]
	}
	
	# Priors for regression coefficients (cP)
	cP <- matrix(NA, nrow=2, ncol=2)
	cP[1,1] <- priors["gppMu"] # prior mean of GPP coefficient (C[1,1]*PAR=GPP)
	cP[1,2] <- priors["gppSig2"] # prior variance of GPP coefficient
	cP[2,1] <- priors["rMu"] # prior mean of R coefficient (C[2,1]*log(Temp)=R)
	cP[2,2] <- priors["rSig2"] # prior variance of R coefficient


	# Put in final format supplied to jags
	data <- list(Y=do.obs, N=length(do.obs), U=U, kP=kP, cP=cP, satO=do.sat, a0=do.obs[1], Zmix=z.mix)
	params <- c("C", "sigmaV", "sigmaW")

	output <- bayesFit(data, params, mf=modfile)
	return(output)
	
}





