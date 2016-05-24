# 2014-09-01 CJS conversion to JAGS
#                - no model name
#                - C(,20) -> T(,20)
#                - dflat() to dnorm(0, 1E-6)
#                - added u2.W.YoYcopy to improve mixing based on Matt S. suggestion
#                - added u2.W.1copy   to improve mixing based on Matt S. suggestion
#                - added u2.H.1copy   to improve mixing based on Matt S. suggestion
#                - fixed monitoring of *H.1 parameters that are only present in hatch.after or later
#                  JAGS won't monitor these variables unless entries from 1:hatch.after are defined
# 2013-09-06 CJS added initializations for n1, m2, u2.W.YoY, u2.W.1, u2.H.1
#                when these are missing (typically when set to bad).
#                Also removed references to WinBugs
# 2011-05-15 CJS limited etaU to 20 or less
# 2011-01-24 SB  added call to run.windows.openbugs and run.windows.winbugs
# 2010-11-25 CJS output to track progress of burnin and post-burnin phases
# 2010-04-26 CJS fixed problem with init.logiP that failed when n1=m2 logit=+infinity and lm() failed
# 2009-12-05 CJS added title to argument list
# 2009-12-01 CJS added openbugs/winbugs directory to argument list

TimeStratPetersenDiagErrorWHSteel <-
  function(title, prefix,
           time, n1, m2,
           u2.W.YoY, u2.W.1, u2.H.1,
           hatch.after=NULL,
           logitP.cov,
           n.chains=3, n.iter=200000, n.burnin=100000, n.sims=2000,
           tauU.alpha=1, tauU.beta=.05, taueU.alpha=1, taueU.beta=.05,
           mu_xiP=logit( sum(m2,na.rm=TRUE)/sum(n1,na.rm=TRUE)),
           tau_xiP=1/var( logit((m2+.5)/(n1+1)),na.rm=TRUE),
           tauP.alpha=.001, tauP.beta=.001,
           debug=FALSE, debug2=FALSE,
	   engine=c('jags','openbugs')[1], 
           InitialSeed){

set.seed(InitialSeed)  # set prior to initial value computations


#
#  Fit the smoothed time-Stratified Petersen estimator with Diagonal recoveries (i.e. no recoveries
#  outside stratum of release), error in the smoothed U curve, and separating wild vs hatchery stocks
#  for steelhead where 100% of hatchery fish are cliped
#
#  This routine assumes that the strata are time (e.g. weeks).
#  In each stratum n1 fish are released (with marks). These are ususally
#     captured fish that are marked, transported upstream, and released.
#     These fish are used only to estimate the recapture rate downstream.
#  Of the n1 fish released, m2 fish are recaptured in the same stratum (e.g. week) of release.
#     There is a related function that allows fish to be recaptured in subsequent weeks.
#
#  There are 3 distinct populations in the stream
#    u2.W.YoY - young of year wild populations
#    u2.W.1   - wild populations of age 1+
#    u2.H.1   - hatchery population of age 1 that appears after hatch.after

#  Input
#      prefix - prefix for file name for initial plot of U's
#      time  - the stratum number
#      n1    - vector of number of fish released in stratum i
#      m2    - vector of number of fish recovered in stratum i (EXCLUDING recaps)
#      u2.W.YoY  - vector of number of wild YoY captured in stratum i
#      u2.W.1    - vector of number of wild YoY captured in stratum i
#      u2.H.1    - vector of number of wild YoY captured in stratum i
#      hatch.after - point AFTER which the hatchery fish are released.
#      logitP.cov - covariates for logit(P)=X beta.logitP


#  This routine makes a call to the MCMC sampler to fit the model and then gets back the
#  coda files for the posteriour distribution.

## Set working directory to current directory (we should allow users to select this)
working.directory <- getwd()

## Define paths for the model, data, and initial value files
model.file <- file.path(working.directory, "model.txt")
data.file <- file.path(working.directory,"data.txt")
init.files <- file.path(working.directory,
                       paste("inits", 1:n.chains,".txt", sep = ""))

# Save the Bugs progam to the model.txt file
#
sink(model.file)  # NOTE: NO " allowed in model as this confuses the cat command
cat("
model {
# Time Stratified Petersen with Diagonal recapture (no spillover in subsequent weeks or marked fish)
#    and allowing for error in the smoothed U curve with separation of wild and hatchery fish
#    for steel head stocks in the Trinity River
# Each of the populations are fit using a SINGLE spline curve as this should be flexible
#    enough to capture the individual behaviours

#  Data input:
#      Nstrata - number of strata
#      n1         - number of marked fish released
#      m2         - number of marked fish recaptured
#      u2.W.YoY   - number of wild     YoY fish captured
#      u2.W.1     - number of wild     1+  fish captured
#      u2.H.1     - number of hatchery 1+  fish captured
#      logitP.cov   - covariates for logitP
#      NlogitP.cov  - number of logitP covariates
#      SplineDesign.W.YoY- wild YoY    fish spline design matrix of size [Nstrata, maxelement of n.b.notflat.W.YoY]
#      SplineDesign.W.1  - wild 1+     fish spline design matrix of size [Nstrata, maxelement of n.b.notflat.W.1  ]
#      SplineDesign.H.1  - hatchery 1+ fish spline design matrix of size [Nstrata, maxelement of n.b.notflat.H.1]
#                   This is set up prior to the call.
#      b.flat.W.YoY   - vector of strata indices where the prior for the b's will be flat for wild YoY fish
#      b.flat.W.1     - vector of strata indices where the prior for the b's will be flat for wild 1+  fish
#      b.flat.H.1     - vector of strata indices where the prior for the b's will be flat for hatchery 1+ fish
#                 this is normally the first two weeks of each spline segment
#      n.b.flat.W.YoY - number of b coefficients that have a flat prior - wild YoY fish
#      n.b.flat.W.1   - number of b coefficients that have a flat prior - wild 1+  fish
#      n.b.flat.H.1   - number of b coefficients that have a flat prior - hatchery fish
#      b.notflat.W.YoY - vector of strata indices where difference in coefficients is modelled - wild YoY fish
#      b.notflat.W.1   - vector of strata indices where difference in coefficients is modelled - wild 1+  fish
#      b.notflat.H.1   - vector of strata indices where difference in coefficients is modelled - hatchery 1+ fish
#      n.b.notflat.W.YoY - number of b coefficients that do not have a flat prior - wild YoY fish
#      n.b.notflat.W.1   - number of b coefficients that do not have a flat prior - wild 1+  fish
#      n.b.notflat.H.1   - number of b coefficients that do not have a flat prior - hatchery 1+ fish
#      tauU.alpha, tauU.beta   - parameters for prior on tauU
#      taueU.alpha, taueU.beta - parameters for prior on taueU
#      mu_xiP, tau_xiP  - parameters for prior on mean logit(P)'s [The intercept term]
#                       - the other beta terms are given a prior of a N(mu=0, variance=1000)
#      tauP.alpha, tauP.beta - parameter for prior on tauP (residual variance of logit(P)'s after adjusting for
#                         covariates)
#
#  Parameters of the model are:
#      p[i]
#       logitP[i]  = logit(p[i]) = logitP.cov*beta.logitP
#         The first beta.logitP has a prior from N(xiP, tauP)
#            and xiP and tauP are currently set internally
#         The remaining beta's are assigned a wider prior N(mu=0, var=1000).
#      U.W.YoY[i] - number of unmarked wild YoY    fish passing stratam i in population
#      U.W.1  [i] - number of unmarked wild 1+     fish passing stratam i in population
#      U.H.1  [i] - number of unmarked hatchery 1+ fish passing stratum i in population
#       etaU.W.YoY[i]  = log(U.W.YoY[i])
#       etaU.W.1  [i]  = log(U.W.1  [i])
#       etaU.H.1  [i]  = log(U.H.1  [i])
#         which comes from spline with parameters bU.*[1... Knots+q] + error term eU.*[i]

   ##### Fit the spline for W.YoY - this covers the entire experiment ######
   for(i in 1:Nstrata){
        logUne.W.YoY[i] <- inprod(SplineDesign.W.YoY[i,1:n.bU.W.YoY],bU.W.YoY[1:n.bU.W.YoY])  # spline design matrix * spline coeff
", fill=TRUE)
sink()  # Temporary end of saving bugs program
if(tolower(engine)=="jags") {
   sink("model.txt", append=TRUE)
   cat("
        etaU.W.YoY[i] ~ dnorm(logUne.W.YoY[i], taueU)T(,20)          # add random error
   ",fill=TRUE)
   sink()
}
if(tolower(engine) %in% c("openbugs")) {
   sink("model.txt", append=TRUE)
   cat("
        etaU.W.YoY[i] ~ dnorm(logUne.W.YoY[i], taueU)C(,20)          # add random error
   ",fill=TRUE)
   sink()
}
   sink("model.txt", append=TRUE)
   cat("
        eU.W.YoY[i] <- etaU.W.YoY[i] - logUne.W.YoY[i]
   }
   ##### Fit the spline for W.1 -   this covers the entire experiment ######
   for(i in 1:Nstrata){
        logUne.W.1[i] <- inprod(SplineDesign.W.1[i,1:n.bU.W.1],bU.W.1[1:n.bU.W.1])  # spline design matrix * spline coeff
   ", fill=TRUE)
sink()  # Temporary end of saving bugs program
if(tolower(engine)=="jags") {
   sink("model.txt", append=TRUE)
   cat("
        etaU.W.1[i] ~ dnorm(logUne.W.1[i], taueU)T(,20)              # add random error
   ",fill=TRUE)
   sink()
}
if(tolower(engine) %in% c("openbugs")) {
   sink("model.txt", append=TRUE)
   cat("
        etaU.W.1[i] ~ dnorm(logUne.W.1[i], taueU)C(,20)              # add random error
   ",fill=TRUE)
   sink()
}
   sink("model.txt", append=TRUE)
   cat("
        eU.W.1[i] <- etaU.W.1[i] - logUne.W.1[i]
   }
   ##### Fit the spline for hatchery fish - these fish only enter AFTER hatch.after ######
   for(i in (hatch.after+1):Nstrata){
        logUne.H.1[i] <- inprod(SplineDesign.H.1[i,1:n.bU.H.1],bU.H.1[1:n.bU.H.1])  # spline design matrix * spline coeff
   ", fill=TRUE)
sink()  # Temporary end of saving bugs program
if(tolower(engine)=="jags") {
   sink("model.txt", append=TRUE)
   cat("
        etaU.H.1[i] ~ dnorm(logUne.H.1[i], taueU)T(,20)              # add random error
   ",fill=TRUE)
   sink()
}
if(tolower(engine) %in% c("openbugs")) {
   sink("model.txt", append=TRUE)
   cat("
        etaU.H.1[i] ~ dnorm(logUne.H.1[i], taueU)C(,20)              # add random error
   ",fill=TRUE)
   sink()
}
   sink("model.txt", append=TRUE)
   cat("
        eU.H.1[i] <- etaU.H.1[i] - logUne.H.1[i]
   }

   ##### Model the capture probabilities #####
   for(i in 1:hatch.after){
        mu.logitP[i] <- inprod(logitP.cov[i,1:NlogitP.cov], beta.logitP[1:NlogitP.cov])
        ## logitP[i] ~ dnorm(mu.logitP[i],tauP)

        # Use the u2.W.YoYcopy to break the cycle (in OpenBugs/Jags) and improve mixing
        mu.epsilon[i] <- mu.logitP[i] - log(u2.W.YoYcopy[i] + u2.W.1copy[i] + 1) +
                    log(exp(etaU.W.YoY[i]) + exp(etaU.W.1[i]))
        epsilon[i] ~ dnorm(mu.epsilon[i],tauP)

        logitP[i] <-  log(u2.W.YoYcopy[i] + u2.W.1copy[i] + 1) -
                    log(exp(etaU.W.YoY[i]) + exp(etaU.W.1[i])) + epsilon[i]
   }
   for(i in (hatch.after+1):Nstrata){
        mu.logitP[i] <- inprod(logitP.cov[i,1:NlogitP.cov], beta.logitP[1:NlogitP.cov])
        ## logitP[i] ~ dnorm(mu.logitP[i],tauP)

        mu.epsilon[i] <- mu.logitP[i] - log(u2.W.YoYcopy[i] + u2.W.1copy[i] + u2.H.1copy[i] + 1) +
                    log(exp(etaU.W.YoY[i]) + exp(etaU.W.1[i]) + exp(etaU.H.1[i]))

        epsilon[i] ~ dnorm(mu.epsilon[i],tauP)

        logitP[i] <-  log(u2.W.YoYcopy[i] + u2.W.1copy[i] + u2.H.1copy[i] + 1) -
                    log(exp(etaU.W.YoY[i]) + exp(etaU.W.1[i]) + exp(etaU.H.1[i])) + epsilon[i]
   }


   ##### Hyperpriors #####
   ## Run size - wild and hatchery fish - flat priors
      ", fill=TRUE)
sink()  # Temporary end of saving bugs program
if(tolower(engine)=="jags") {
   sink("model.txt", append=TRUE)
   cat("
   for(i in 1:n.b.flat.W.YoY){
      bU.W.YoY[b.flat.W.YoY[i]] ~ dnorm(0, 1E-6)
   }
   for(i in 1:n.b.flat.W.1){
      bU.W.1[b.flat.W.1[i]] ~ dnorm(0, 1E-6)
   }
   for(i in 1:n.b.flat.H.1){
      bU.H.1[b.flat.H.1[i]] ~ dnorm(0, 1E-6)
   }
   ",fill=TRUE)
   sink()
}
if(tolower(engine) %in% c("openbugs")) {
   sink("model.txt", append=TRUE)
   cat("
   for(i in 1:n.b.flat.W.YoY){
      bU.W.YoY[b.flat.W.YoY[i]] ~ dflat()
   }
   for(i in 1:n.b.flat.W.1){
      bU.W.1[b.flat.W.1[i]] ~ dflat()
   }
   for(i in 1:n.b.flat.H.1){
      bU.H.1[b.flat.H.1[i]] ~ dflat()
   }
   ",fill=TRUE)
   sink()
}
   sink("model.txt", append=TRUE)
   cat("
   ## Run size - priors on the difference for wild and hatchery fish
   for(i in 1:n.b.notflat.W.YoY){
      xiU.W.YoY[b.notflat.W.YoY[i]] <- 2*bU.W.YoY[b.notflat.W.YoY[i]-1] - bU.W.YoY[b.notflat.W.YoY[i]-2]
      bU.W.YoY [b.notflat.W.YoY[i]] ~ dnorm(xiU.W.YoY[b.notflat.W.YoY[i]],tauU)
   }
   for(i in 1:n.b.notflat.W.1){
      xiU.W.1[b.notflat.W.1[i]] <- 2*bU.W.1[b.notflat.W.1[i]-1] - bU.W.1[b.notflat.W.1[i]-2]
      bU.W.1 [b.notflat.W.1[i]] ~ dnorm(xiU.W.1[b.notflat.W.1[i]],tauU)
   }
   for(i in 1:n.b.notflat.H.1){
      xiU.H.1[b.notflat.H.1[i]] <- 2*bU.H.1[b.notflat.H.1[i]-1] - bU.H.1[b.notflat.H.1[i]-2]
      bU.H.1 [b.notflat.H.1[i]] ~ dnorm(xiU.H.1[b.notflat.H.1[i]],tauU)
   }

   tauU ~ dgamma(tauU.alpha,tauU.beta)  # Notice reduction from .0005 (in thesis) to .05
   sigmaU <- 1/sqrt(tauU)
   taueU ~ dgamma(taueU.alpha,taueU.beta) # dgamma(100,.05) # Notice reduction from .0005 (in thesis) to .05
   sigmaeU <- 1/sqrt(taueU)

   ## Capture probabilities. The logit(p[i]) are n(logitP.cov*beta.logitP.cov, sigmaP**2)
   beta.logitP[1] ~ dnorm(mu_xiP,tau_xiP) # first term is usually an intercept
   for(i in 2:NlogitP.cov){
      beta.logitP[i] ~ dnorm(0, .01)   # rest of beta terms are normal 0 and a large variance
   }
   beta.logitP[NlogitP.cov+1] ~ dnorm(0, .01) # dummy so that covariates of length 1 function properly
   tauP ~ dgamma(tauP.alpha,tauP.beta)
   sigmaP <- 1/sqrt(tauP)

   ##### Likelihood contributions #####
   ## Number of marked fish recovered ##
   for(i in 1:Nstrata){
      logit(p[i]) <- logitP[i]       # convert from logit scale
      m2[i] ~ dbin(p[i],n1[i])     # recovery of marked fish
   }

   ## captures of wild YoY and wild 1+  fish
   for(i in 1:Nstrata){
      U.W.YoY[i] <- round(exp(etaU.W.YoY[i]))       # convert from log scale
      U.W.1  [i] <- round(exp(etaU.W.1  [i]))       # convert from log scale
      u2.W.YoY[i] ~ dbin(p[i],U.W.YoY[i])
      u2.W.1  [i] ~ dbin(p[i],U.W.1  [i])
   }

   ## captures of hatchery fish - these can only occur AFTER hatch.after
   for(i in (hatch.after+1):Nstrata){
      U.H.1[i] <- round(exp(etaU.H.1[i]))      # convert from log scale
      u2.H.1[i] ~ dbin(p[i], U.H.1[i])
   }

   ##### Derived Parameters #####
   Utot.W.YoY <- sum( U.W.YoY[1:Nstrata])              # Total number of unmarked fish - wild YoY
   Utot.W.1   <- sum( U.W.1  [1:Nstrata])              # Total number of unmarked fish - wild 1+
   Utot.H.1   <- sum( U.H.1  [(hatch.after+1):Nstrata])# Total number of unmarked fish - hatchery 1+
   Utot   <- Utot.W.YoY + Utot.W.1 + Utot.H.1                   # Grand total number of fish

   # Because JAGES does not properly monitory partially defined vectors (see Section 2.5 of the JAGES user manual)
   # we need to add dummy distribution for the parameters of interest prior to the hatchery fish arriving.
   # This is not needed in OpenBugs who returns the subset actually monitored, but we add this to be consistent
   # among the two programs
   for(i in 1:hatch.after){
      U.H.1[i]      ~ dnorm(0,1)  # These are complete arbitrary and never gets updated
      etaU.H.1[i]   ~ dnorm(0,1)
      logUne.H.1[i] ~ dnorm(0,1)
      eU.H.1[i]     ~ dnorm(0,1)
   }
}  # end of model

", fill=TRUE)
sink()  # End of saving the Bugs program

Nstrata <- length(n1)

# Make a copy of u2.W.YoY/u2.W.1/u2.H.1 to improve mixing in the MCMC model
u2.W.YoYcopy <- spline(x=1:Nstrata, y=u2.W.YoY, xout=1:Nstrata)$y
u2.W.YoYcopy <- round(u2.W.YoYcopy) # round to integers
u2.W.1copy   <- spline(x=1:Nstrata, y=u2.W.1,   xout=1:Nstrata)$y
u2.W.1copy   <- round(u2.W.1copy) # round to integers

# similarly make a copy of u2.H.1 to improve mixing in the MCMC model
# notice that hatchery fish occur at hatch.after or later
u2.H.1copy <- u2.H.1 * 0
u2.H.1copy[hatch.after:Nstrata] <- spline(x=hatch.after:Nstrata, y=u2.H.1[hatch.after:Nstrata], xout=hatch.after:Nstrata)$y
u2.H.1copy <- round(u2.H.1copy) # round to integers




datalist <- list("Nstrata", "n1", "m2", 
		 "u2.W.YoY", "u2.W.YoYcopy",
		 "u2.W.1",   "u2.W.1copy",
		 "u2.H.1",   "u2.H.1copy",
		 "hatch.after",
                 "logitP.cov", "NlogitP.cov",
                 "SplineDesign.W.YoY",
                 "b.flat.W.YoY", "n.b.flat.W.YoY", "b.notflat.W.YoY", "n.b.notflat.W.YoY", "n.bU.W.YoY",
                 "SplineDesign.W.1",
                 "b.flat.W.1",   "n.b.flat.W.1",   "b.notflat.W.1",   "n.b.notflat.W.1",   "n.bU.W.1",
                 "SplineDesign.H.1",
                 "b.flat.H.1", "n.b.flat.H.1", "b.notflat.H.1", "n.b.notflat.H.1", "n.bU.H.1",
                 "tauU.alpha", "tauU.beta", "taueU.alpha", "taueU.beta",
                 "mu_xiP", "tau_xiP", "tauP.alpha", "tauP.beta")

parameters <- c("logitP", "beta.logitP", "tauP", "sigmaP",
                "bU.W.YoY", "bU.W.1", "bU.H.1", "tauU", "sigmaU",
                "eU.W.YoY", "eU.W.1", "eU.H.1", "taueU", "sigmaeU",
                "Utot.W.YoY", "Utot.W.1", "Utot.H.1", "Utot", "logUne.W.YoY", "logUne.W.1", "logUne.H.1",
                "etaU.W.YoY", "etaU.W.1", "etaU.H.1", "U.W.YoY", "U.W.1", "U.H.1")
if( any(is.na(m2))) {parameters <- c(parameters,"m2")} # monitor in case some bad data where missing values present
if( any(is.na(u2.W.YoY))) {parameters <- c(parameters,"u2.W.YoY")}
if( any(is.na(u2.W.1)))   {parameters <- c(parameters,"u2.W.1")}
if( any(is.na(u2.H.1)))   {parameters <- c(parameters,"u2.H.1")}


# Now to create the initial values, and the data prior to call to the MCMC sampler

Nstrata <- length(n1)

# Estimate number of wild and hatchery fish based on clip rate
u2.W.YoY[is.na(u2.W.YoY)] <- 1  # in case of missing values
u2.W.1  [is.na(u2.W.1)]   <- 1  # in case of missing values
u2.H.1  [is.na(u2.H.1)]   <- 1  # in case of missing values

avg.P <- sum(m2,na.rm=TRUE)/sum(n1, na.rM=TRUE)
Uguess.W.YoY <- pmax((u2.W.YoY+1)*(n1+2)/(m2+1), u2.W.YoY/avg.P, 1,na.rm=TRUE)  # try and keep Uguess larger than observed values
Uguess.W.1   <- pmax((u2.W.1  +1)*(n1+2)/(m2+1), u2.W.1  /avg.P, 1, na.rm=TRUE)  # try and keep Uguess larger than observed values
Uguess.H.1   <- pmax((u2.H.1  +1)*(n1+2)/(m2+1), u2.H.1  /avg.P, 1, na.rm=TRUE)  # try and keep Uguess larger than observed values
Uguess.H.1[1:hatch.after] <- 0   # no hatchery fish prior to release from hatchery

# create the B-spline design matrix for wild and hatchery fish
# The design matrix for hatchery fill will still have rows corresponding to entries PRIOR to
#    the hatchery release but these are never used in the winbugs fitting routines
# There is a separate (single) spline for hatchery and wild fish with NO breakpoints
# The first two coefficient have a flat prior and the rest of the coefficients are modelled using
#    differences between the succesive coefficients

# Wild YoY  fish. This covers the entire experiment.
SplineDegree <- 3           # Degree of spline between occasions
knots <- seq(4,Nstrata,4)/(Nstrata+1) # a knot roughly every 4th stratum
SplineDesign.W.YoY <- bs((1:Nstrata)/(Nstrata+1), knots=knots, degree=SplineDegree, intercept=TRUE, Boundary.knots=c(0,1))
b.flat.W.YoY      <- c(1,2)
b.notflat.W.YoY   <- 3:(ncol(SplineDesign.W.YoY))
n.b.flat.W.YoY    <- length(b.flat.W.YoY)
n.b.notflat.W.YoY <- length(b.notflat.W.YoY)
n.bU.W.YoY        <- n.b.flat.W.YoY + n.b.notflat.W.YoY
init.bU.W.YoY    <- lm(log(Uguess.W.YoY+1) ~ SplineDesign.W.YoY-1)$coefficients  # initial values for spline coefficients

# Wild 1+  fish. This covers the entire experiment.
SplineDegree <- 3           # Degree of spline between occasions
knots <- seq(4,Nstrata,4)/(Nstrata+1) # a knot roughly every 4th stratum
SplineDesign.W.1   <- bs((1:Nstrata)/(Nstrata+1), knots=knots, degree=SplineDegree, intercept=TRUE, Boundary.knots=c(0,1))
b.flat.W.1        <- c(1,2)
b.notflat.W.1     <- 3:(ncol(SplineDesign.W.1))
n.b.flat.W.1      <- length(b.flat.W.1)
n.b.notflat.W.1   <- length(b.notflat.W.1)
n.bU.W.1          <- n.b.flat.W.1 + n.b.notflat.W.1
init.bU.W.1       <- lm(log(Uguess.W.1+1) ~ SplineDesign.W.1-1)$coefficients  # initial values for spline coefficients

# hatchery fish. Notice they can only enter AFTER hatch.after, The spline design matrix still has rows
# of zero for 1:hatch.after to make it easier in Bugs
SplineDegree <- 3           # Degree of spline between occasions
knots <- (seq((hatch.after+4),Nstrata-1,4)-hatch.after)/(Nstrata-hatch.after+1) # a knot roughly every 4th stratum
SplineDesign.H.1 <- bs((1:(Nstrata-hatch.after))/(Nstrata-hatch.after+1), knots=knots, degree=SplineDegree, intercept=TRUE, Boundary.knots=c(0,1))
b.flat.H.1      <- c(1,2)
b.notflat.H.1   <- 3:(ncol(SplineDesign.H.1))
n.b.flat.H.1    <- length(b.flat.H.1)
n.b.notflat.H.1 <- length(b.notflat.H.1)
n.bU.H.1        <- n.b.flat.H.1 + n.b.notflat.H.1
init.bU.H.1     <- lm(log(Uguess.H.1[(hatch.after+1):Nstrata]+1) ~ SplineDesign.H.1-1)$coefficients  # initial values for spline coefficients
# patch up the initial rows of the spline design matrix
SplineDesign.H.1 <- rbind(matrix(0,nrow=hatch.after, ncol=ncol(SplineDesign.H.1)), SplineDesign.H.1)


#browser()

# get the logitP=logit(P) covariate matrix ready
logitP.cov <- as.matrix(logitP.cov)
NlogitP.cov <- ncol(as.matrix(logitP.cov))


#  initial values for the parameters of the model

## init.vals <- function(){
##    # Initial values for the probability of capture
##    init.logitP <- pmax(-10,pmin(10,logit((m2+1)/(n1+2))))  # initial capture rates based on observed recaptures
##    init.logitP[is.na(init.logitP)] <- -2         # those cases where initial probability is unknown
##    init.beta.logitP <- as.vector(lm( init.logitP ~ logitP.cov-1)$coefficients)
##    init.beta.logitP[init.beta.logitP=NA] <- 0
##    init.beta.logitP <- c(init.beta.logitP, 0)   # add one extra element so that single beta is still written as a
##                                              # vector in the init files etc.
##    init.tauP <- 1/var(init.logitP, na.rm=TRUE)     # 1/variance of logit(p)'s (ignoring the covariates for now)

##    # inital values for the spline coefficients
##    init.bU.W.YoY   <- lm(log(Uguess.W.YoY+1) ~ SplineDesign.W.YoY-1)$coefficients  # initial values for spline coefficients
##    init.bU.W.1     <- lm(log(Uguess.W.1  +1) ~ SplineDesign.W.1  -1)$coefficients  # initial values for spline coefficients
##    init.bU.H.1     <- lm(log(Uguess.H.1[(hatch.after+1):Nstrata]+1) ~ SplineDesign.H.1[(hatch.after+1):Nstrata,]-1)$coefficients  # initial values for spline coefficients

##    init.eU.W.YoY   <- as.vector(log(Uguess.W.YoY+1)-SplineDesign.W.YoY%*%init.bU.W.YoY)  # error terms set as differ between obs and pred
##    init.eU.W.1     <- as.vector(log(Uguess.W.1  +1)-SplineDesign.W.1  %*%init.bU.W.1)  # error terms set as differ between obs and pred
##    init.eU.H.1     <- as.vector(log(Uguess.H.1  +1)-SplineDesign.H.1  %*%init.bU.H.1)  # error terms set as differ between obs and pred

##    init.etaU.W.YoY <- log(Uguess.W.YoY+1)
##    init.etaU.W.1   <- log(Uguess.W.1  +1)
##    init.etaU.H.1   <- log(Uguess.H.1  +1)
##    init.etaU.H.1[1:hatch.after] <- NA  # these are never used.

##    # variance of spline difference (use only the wild fish to initialize)
##    sigmaU <- sd( init.bU.W.YoY[b.notflat.W.YoY]-2*init.bU.W.YoY[b.notflat.W.YoY-1]+init.bU.W.YoY[b.notflat.W.YoY-2], na.rm=TRUE)
##    init.tauU <- 1/sigmaU^2

##    # variance of error in the U' over and above the spline fit (use only the wild fish to initialize)
##    sigmaeU <- sd(init.eU.W.YoY, na.rm=TRUE)
##    init.taueU <- 1/sigmaeU^2

##    # initialize the u2.* where missing
##    init.u2.W.YoY    <- u2.W.YoY
##    init.u2.W.YoY[ is.na(u2.W.YoY)] <- 100
##    init.u2.W.YoY[!is.na(u2.W.YoY)] <- NA

##    init.u2.W.1      <- u2.W.1
##    init.u2.W.1  [ is.na(u2.W.1)] <- 100
##    init.u2.W.1  [!is.na(u2.W.1)] <- NA

##    init.u2.H.1      <- u2.H.1
##    init.u2.H.1  [ is.na(u2.H.1)] <- 100
##    init.u2.H.1  [!is.na(u2.H.1)] <- NA

##    list(logitP=init.logitP, beta.logitP=init.beta.logitP, tauP=init.tauP,
##         bU.W.YoY=init.bU.W.YoY,
##         bU.W.1  =init.bU.W.1,
##         bU.H.1  =init.bU.H.1,
##         tauU=init.tauU, taueU=init.taueU,
##         etaU.W.YoY=init.etaU.W.YoY, etaU.W.1=init.etaU.W.1, etaU.H.1=init.etaU.H.1)
## }
## #browser()

##  initial values for the parameters of the model
init.vals <- genInitVals(model="TSPDE-WHsteel",
                         n1=n1,
                         m2=m2,
                         u2=list(W.YoY=u2.W.YoY,W.1=u2.W.1,H.1=u2.H.1),
                         logitP.cov=logitP.cov,
                         SplineDesign=list(W.YoY=SplineDesign.W.YoY,W.1=SplineDesign.W.1,H.1=SplineDesign.H.1),
                         hatch.after=hatch.after,
                         n.chains=n.chains)

## Generate data list
data.list <- list()
for(i in 1:length(datalist)){
  data.list[[length(data.list)+1]] <-get(datalist[[i]])
}
names(data.list) <- as.list(datalist)

# Call the MCMC sampler
results <- run.MCMC(modelFile=model.file,
                        dataFile=data.file,
                        dataList=data.list,
                        initFiles=init.files,
                        initVals=init.vals,
                        parameters=parameters,
                        nChains=n.chains,
                        nIter=n.iter,
                        nBurnin=n.burnin,
                        nSims=n.sims,
                        overRelax=FALSE,
                        initialSeed=InitialSeed,
                        working.directory=working.directory,
			engine=engine,
                        debug=debug)

return(results)
}
