# 2013-12-31 CJS conversion to JAGS
#                - no model name
#                - C(,20) -> T(,20)
#                - dflat() to dnorm(0, 1E-6)
#                - added u2.N.1copy to improve mixing based on Matt S. suggestion   ????? - not done here?
#                - added u2.A.1copy to improve mixing based on Matt S. suggestion
#                - added u2.N.YoYcopy to improve mixing based on Matt S. suggestion
#                - added u2.A.YoYcopy to improve mixing based on Matt S. suggestion
#                - fixed monitoring of *H parameters that are only present in hatch.after or later
#                  JAGS won't monitor these variables unless entries from 1:hatch.after are defined
# 2011-05-15 CJS limited etaU to 20 or less
# 2011-01-24 SB  added call to run.windows.openbugs and run.windows.winbugs
# 2010-11-25 CJS add output to track progress of sampling through burnin and post-burnin
# 2010-04-26 CJS fixed problem where init.logitP failed when n1=m2 (logit=infinite) and lm() failed.
# 2010-03-29 CJS Created first release

# This DIFFERS from the TimeStratPetersenDiagErrorWHChinook routine in the following ways.
#   YoY chinook are separated from age 1 chinook
#   The wild YoY chinook are present in the stream with NO AD clips until the hatchery fish arrive 
#   The Age 1 chinook (from last year) are still present for the entire experiment with some of them
#   having ad-fin clips using the clip-rate from last year.
# The n1/m2 recapture portion is assumed to be common to both ages of fish (a doubtful assumption?)

TimeStratPetersenDiagErrorWHChinook2 <-
       function(title, prefix, time, n1, m2, 
                u2.A.YoY, u2.N.YoY, u2.A.1, u2.N.1,
                hatch.after.YoY=NULL, 
                clip.frac.H.YoY=.25, clip.frac.H.1 = .25,
                logitP.cov,
                n.chains=3, n.iter=200000, n.burnin=100000, n.sims=2000, 
                tauU.alpha=1, tauU.beta=.05, taueU.alpha=1, taueU.beta=.05,
                mu_xiP=logit(sum(m2,na.rm=TRUE)/sum(n1,na.rm=TRUE)), 
                tau_xiP=1/var(logit((m2+.5)/(n1+1)), na.rm=TRUE), 
                tauP.alpha=.001, tauP.beta=.001, 
                debug=FALSE, debug2=FALSE,
		engine=c('jags','openbugs')[1],
                InitialSeed){

set.seed(InitialSeed)  # set prior to initial value computations


#
#  Fit the smoothed time-Stratified Petersen estimator with Diagonal recoveries (i.e. no recoveries
#  outside stratum of release), error in the smoothed U curve, and separating wild vs hatchery stocks
#  for chinook. 
#
#  This routine assumes that the strata are time (e.g. weeks).
#  In each stratum n1 fish are released (with marks). These are ususally
#     captured fish that are marked, transported upstream, and released.
#     These fish are used only to estimate the recapture rate downstream.
#  Of the n1 fish released, m2 fish are recaptured in the same stratum (e.g. week) of release.
#     There is a related function that allows fish to be recaptured in subsequent weeks.
#
#  Both YoY and age 1 fish are present in the stream.
#  The traps capture u2.A.YoY (YoY with adipose fin clipped) and u2.N.YoY (YoY not clipped)
#     are newly captured in stratum i. 
#     Prior to the hatch.after.YoY, there are no ad-fin clipped YoY fish and all YoY fish captured
#     are assumed to be wild. The clip-fraction of the hatchery fish is clip.frac.H.YoY
#  The traps also capture u2.A.1 (age 1 with ad-fin clipped) and u2.N.1 (age 1 with no ad-fin clipped)
#     which represent fish from LAST year that have residualized in the stream. These are a mixture of
#     wild and hatchery fish. The clip rate for these fish from last year is clip.frac.H.1
#
#  All wild fish are NOT ad-clipped.
#  Only a fraction of hatchery fish are ad-clipped. It is assumed that the fraction of ad-clipped
#     hatchery fish is constant over the life of the run.       
# 
#  Input
#      prefix - prefix for file name for initial plot of U's
#      time  - the stratum number
#      n1    - vector of number of fish released in stratum i
#      m2    - vector of number of fish recovered in stratum i (EXCLUDING recaps)
#      u2.A.YoY  - vector of number of ad-clipped  YoY unmarked fish captured in stratum i
#      u2.N.YoY  - vector of number of non-clipped YoY unmarked fish captured in stratum i
#      u2.A.1    - vector of number of ad-clipped  Age1 unmarked fish captured in stratum i
#      u2.N.1    - vector of number of non-clipped Age1 unmarked fish captured in stratum i
#      hatch.after.YoY - point AFTER which the YoY hatchery fish are released.
#      clip.frac.H.YoY - what fraction of hatchery fish are clipped at YoY
#      clip.frac.H.1   - what fraction of hatchery fish are clipped as YoY LAST YEAR who are now age 1
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
# Each of the wild and hatchery populations are fit using a SINGLE spline curve as this should be flexible 
#    enough to capture the individual behaviours

#  Data input:
#      Nstrata - number of strata
#      n1         - number of marked fish released
#      m2         - number of marked fish recaptured
#      u2.A.YoY   - number of adclipped   YoY unmarked fish captured (must be hatchery fish).
#      u2.N.YoY   - number of non-clipped YoY unmarked fish captured (wild + hatchery fish)
#      u2.A.1     - number of adclipped  Age1 unmarked fish captured (must be hatchery fish from last year)
#      u2.A.1     - number of adclipped  Age1 unmarked fish captured (wild + hatchery fish  from last year)
#      clip.frac.H.YoY- what fraction of YoY hatchery fish are clipped 
#      clip.frac.H.1  - what fraction of Age1 hatcery fish are clipped (from last years hatchery release)
#      logitP.cov   - covariates for logitP
#      NlogitP.cov  - number of logitP covariates
#      SplineDesign.W.YoY- YoY wildfish spline design matrix of size [Nstrata, maxelement of n.b.notflat.W.YoY]
#      SplineDesign.H.YoY- YoY hatchery spline design matrix of size [Nstrata, maxelement of n.b.notflat.H.YoY]
#      SplineDesign.W.1  - Age1 wildfish spline design matrix of size [Nstrata, maxelement of n.b.notflat.W.1]
#      SplineDesign.H.1  - Age1 hatchery spline design matrix of size [Nstrata, maxelement of n.b.notflat.H.1]
#                   These design matrices are set up prior to the call.
#      b.flat.W.YoY   - vector of strata indices where the prior for the b's will be flat for YoY wild     fish
#      b.flat.H.YoY   - vector of strata indices where the prior for the b's will be flat for YoY hatchery fish
#      b.flat.W.1     - vector of strata indices where the prior for the b's will be flat for Age1 wild     fish
#      b.flat.H.1     - vector of strata indices where the prior for the b's will be flat for Age1 hatchery fish
#                 this are normally the first two strata of each spline segment
#      n.b.flat.W.YoY - number of b coefficients that have a flat prior - YoY wild     fish
#      n.b.flat.H.YoY - number of b coefficients that have a flat prior - YoY hatchery fish
#      n.b.flat.W.1   - number of b coefficients that have a flat prior - Age1 wild     fish
#      n.b.flat.H.1   - number of b coefficients that have a flat prior - Age1 hatchery fish
#      b.notflat.W.YoY - vector of strata indices where difference in coefficients is modelled - YoY wild     fish
#      b.notflat.H.YoY - vector of strata indices where difference in coefficients is modelled - YoY hatchery fish
#      b.notflat.W.1   - vector of strata indices where difference in coefficients is modelled - Age1 wild     fish
#      b.notflat.H.1   - vector of strata indices where difference in coefficients is modelled - Age1 hatchery fish
#      n.b.notflat.W.YoY - number of b coefficients that do not have a flat prior - YoY wild     fish
#      n.b.notflat.H.YoY - number of b coefficients that do not have a flat prior - YoY hatchery fish
#      n.b.notflat.W.1   - number of b coefficients that do not have a flat prior - Age1 wild     fish
#      n.b.notflat.H.1   - number of b coefficients that do not have a flat prior - Age1 hatchery fish
#      tauU.alpha, tauU.beta   - parameters for prior on tauU
#      taueU.alpha, taueU.beta - parameters for prior on taueU
#      mu_xiP, tau_xiP  - parameters for prior on mean logit(P)'s [The intercept term]
#                       - the other beta terms are given a prior of a N(mu=0, variance=1000)
#      tauP.alpha, tauP.beta - parameter for prior on tauP (residual variance of logit(P)'s after adjusting for
#                         covariates)
#      clip.frac.H.YoY   - what fraction of YoY  hatchery fish are clipped (KNOWN in advance)
#      clip.frac.H.1     - what fraction of Age1 hatchery fish are clipped (KNOWN in advance from last year's releases)
#
#  Parameters of the model are:
#      p[i]
#       logitP[i]  = logit(p[i]) = logitP.cov*beta.logitP
#         The first beta.logitP has a prior from N(xiP, tauP)
#            and xiP and tauP are currently set internally
#         The remaining beta's are assigned a wider prior N(mu=0, var=1000).
#      U.W.YoY[i] - number of YoY unmarked wild     fish passing stratam i in population
#      U.H.YoY[i] - number of YoY unmarked hatchery fish passing stratum i in population
#       etaU.W.YoY[i]  = log(U.W.YoY[i])
#       etaU.H.YoY[i]  = log(U.H.YoY[i])
#         which comes from spline with parameters bU.W.YoY[1... Knots+q] or bU.H.YoY[1... knots+q]
#         + error term eU.W.YoY[i] or eu.H.YoY[i]        
#      U.W.1[i] - number of Age1 unmarked wild     fish passing stratam i in population
#      U.H.1[i] - number of Age1 unmarked hatchery fish passing stratum i in population
#       etaU.W.1[i]  = log(U.W.1[i])
#       etaU.H.1[i]  = log(U.H.1[i])
#         which comes from spline with parameters bU.W.1[1... Knots+q] or bU.H.1[1... knots+q]
#         + error term eU.W.1[i] or eu.H.1[i]        

   ##### Fit the spline for YoY wildfish - this covers the entire experiment ######
   for(i in 1:Nstrata){
        logUne.W.YoY[i] <- inprod(SplineDesign.W.YoY[i,1:n.bU.W.YoY],bU.W.YoY[1:n.bU.W.YoY])  # spline design matrix * spline coeff 
", fill=TRUE)
sink()  # Temporary end of saving bugs program
if(tolower(engine)=="jags") {
   sink("model.txt", append=TRUE)
   cat("
        etaU.W.YoY[i] ~ dnorm(logUne.W.YoY[i], taueU)T(,20)              # add random error
   ",fill=TRUE)
   sink()
}
if(tolower(engine) %in% c("openbugs")) {
   sink("model.txt", append=TRUE)
   cat("
        etaU.W.YoY[i] ~ dnorm(logUne.W.YoY[i], taueU)C(,20)              # add random error
   ",fill=TRUE)
   sink()
}
   sink("model.txt", append=TRUE)
   cat("
        eU.W.YoY  [i] <- etaU.W.YoY[i] - logUne.W.YoY[i]
   }
   ##### Fit the spline for YoY hatchery fish - these fish only enter AFTER hatch.after.YoY ######
   for(i in (hatch.after.YoY+1):Nstrata){
        logUne.H.YoY[i] <- inprod(SplineDesign.H.YoY[i,1:n.bU.H.YoY],bU.H.YoY[1:n.bU.H.YoY])  # spline design matrix * spline coeff 
  ", fill=TRUE)
sink()  # Temporary end of saving bugs program
if(tolower(engine)=="jags") {
   sink("model.txt", append=TRUE)
   cat("
        etaU.H.YoY[i] ~ dnorm(logUne.H.YoY[i], taueU)T(,20)              # add random error
   ",fill=TRUE)
   sink()
}
if(tolower(engine) %in% c("openbugs")) {
   sink("model.txt", append=TRUE)
   cat("
        etaU.H.YoY[i] ~ dnorm(logUne.H.YoY[i], taueU)C(,20)              # add random error
   ",fill=TRUE)
   sink()
}
   sink("model.txt", append=TRUE)
   cat("
        eU.H.YoY  [i] <- etaU.H.YoY[i] - logUne.H.YoY[i]
   }
   ##### Fit the spline for Age1 wildfish - this covers the entire experiment ######
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
        eU.W.1  [i] <- etaU.W.1[i] - logUne.W.1[i]
   }
   ##### Fit the spline for Age1 hatchery fish - this covers the entire experiment because the have residualized from last year
   for(i in 1:Nstrata){
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
        eU.H.1  [i] <- etaU.H.1[i] - logUne.H.1[i]
   }



   ##### Model the capture probabilities #####
   for(i in 1:Nstrata){
        mu.logitP[i] <- inprod(logitP.cov[i,1:NlogitP.cov], beta.logitP[1:NlogitP.cov]) 
        logitP[i] ~ dnorm(mu.logitP[i],tauP)
   }

   ", fill=TRUE)
sink()  # Temporary end of saving bugs program
if(tolower(engine)=="jags") {
   sink("model.txt", append=TRUE)
   cat("
   ##### Hyperpriors #####
   ## Run size - wild and hatchery fish - flat priors
   for(i in 1:n.b.flat.W.YoY){
      bU.W.YoY[b.flat.W.YoY[i]] ~ dnorm(0, 1E-6)
   }
   for(i in 1:n.b.flat.H.YoY){
      bU.H.YoY[b.flat.H.YoY[i]] ~ dnorm(0, 1E-6)
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
   ##### Hyperpriors #####
   ## Run size - wild and hatchery fish - flat priors
   for(i in 1:n.b.flat.W.YoY){
      bU.W.YoY[b.flat.W.YoY[i]] ~ dflat()
   }
   for(i in 1:n.b.flat.H.YoY){
      bU.H.YoY[b.flat.H.YoY[i]] ~ dflat()
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

   ## Run size - priors on the difference for YoY wild and hatchery fish
   for(i in 1:n.b.notflat.W.YoY){
      xiU.W.YoY[b.notflat.W.YoY[i]] <- 2*bU.W.YoY[b.notflat.W.YoY[i]-1] - bU.W.YoY[b.notflat.W.YoY[i]-2]
      bU.W.YoY [b.notflat.W.YoY[i]] ~ dnorm(xiU.W.YoY[b.notflat.W.YoY[i]],tauU)
   }
   for(i in 1:n.b.notflat.H.YoY){
      xiU.H.YoY[b.notflat.H.YoY[i]] <- 2*bU.H.YoY[b.notflat.H.YoY[i]-1] - bU.H.YoY[b.notflat.H.YoY[i]-2]
      bU.H.YoY [b.notflat.H.YoY[i]] ~ dnorm(xiU.H.YoY[b.notflat.H.YoY[i]],tauU)
   }
   ## Run size - priors on the difference for AGE1 wild and hatchery fish
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

   ## captures of YoY wild (unclipped fish) - these are the only fish available upto (and including) hatch.after
   for(i in 1:hatch.after.YoY){
      U.W.YoY[i] <- round(exp(etaU.W.YoY[i]))       # convert from log scale
      u2.N.YoY[i] ~ dbin(p[i],U.W.YoY[i])
   }

   ## captures of YoY hatchery (clipped fish) - these can only occur AFTER hatch.after
   for(i in (hatch.after.YoY+1):Nstrata){
      U.W.YoY[i] <- round(exp(etaU.W.YoY[i]))      # convert from log scale
      U.H.YoY[i] <- round(exp(etaU.H.YoY[i]))      # convert from log scale
      U.clip.YoY[i] ~ dbin(clip.frac.H.YoY, U.H.YoY [i])
      p.temp.YoY[i] <- p[i]*clip.frac.H.YoY
      u2.A.YoY[i] ~ dbin(p.temp.YoY[i], U.H.YoY[i]) # must be hatchery and clipped
   }

   ## captures of YoY wild+hatchery unclipped fish - these can only occur AFTER hatch.after
   for(i in (hatch.after.YoY+1):Nstrata){
      U.noclip.YoY[i] <- U.W.YoY[i] + U.H.YoY[i] - U.clip.YoY[i]
      u2.N.YoY[i] ~ dbin(p[i], U.noclip.YoY[i])
   } 
 
   ## captures of Age1 wild+hatchery  (clipped fish)
   for(i in 1:Nstrata){
      U.W.1[i] <- round(exp(etaU.W.1[i]))      # convert from log scale
      U.H.1[i] <- round(exp(etaU.H.1[i]))      # convert from log scale
      U.clip.1[i] ~ dbin(clip.frac.H.1, U.H.1[i])
      #u2.A[i] ~ dbin(p[i], U.clip.YoY[i])
      p.temp.1[i] <- p[i]*clip.frac.H.1
      u2.A.1[i] ~ dbin(p.temp.1[i], U.H.1[i]) # must be hatchery and clipped
   }

   ## captures of Age1 wild+hatchery unclipped fish
   for(i in 1:Nstrata){
      U.noclip.1[i] <- U.W.1[i] + U.H.1[i] - U.clip.1[i]
      u2.N.1[i] ~ dbin(p[i], U.noclip.1[i])
   } 



   ##### Derived Parameters #####
   Utot.W.YoY <- sum( U.W.YoY[1:Nstrata])              # Total number of YoY unmarked fish - wild
   Utot.H.YoY <- sum( U.H.YoY[(hatch.after.YoY+1):Nstrata])# Total number of YoY unmarked fish - hatchery
   Utot.W.1   <- sum( U.W.1[1:Nstrata])                # Total number of Age1 unmarked fish - wild
   Utot.H.1   <- sum( U.H.1[1:Nstrata])                # Total number of Age1 unmarked fish - hatchery
   Utot.YoY   <- Utot.W.YoY + Utot.H.YoY               # Grand total number of YoY fish
   Utot.1     <- Utot.W.1   + Utot.H.1                 # Grand total number of Age1 fish
   Utot       <- Utot.YoY + Utot.1

   # Because JAGES does not properly monitory partially defined vectors (see Section 2.5 of the JAGES user manual)
   # we need to add dummy distribution for the parameters of interest prior to the hatchery fish arriving.
   # This is not needed in OpenBugs who returns the subset actually monitored, but we add this to be consistent
   # among the two programs
   for(i in 1:hatch.after.YoY){
      U.H.YoY[i]      ~ dnorm(0,1)  # These are complete arbitrary and never gets updated
      etaU.H.YoY[i]   ~ dnorm(0,1)
      logUne.H.YoY[i] ~ dnorm(0,1)
      eU.H.YoY[i]     ~ dnorm(0,1)
   }

}  # end of model

", fill=TRUE)
sink()  # End of saving the Bugs program


datalist <- list("Nstrata", "n1", "m2",
                 "u2.A.YoY", "u2.N.YoY", "u2.A.1", "u2.N.1",
                 "hatch.after.YoY", "clip.frac.H.YoY", "clip.frac.H.1",
                 "logitP.cov", "NlogitP.cov",
                 "SplineDesign.W.YoY",
                 "b.flat.W.YoY", "n.b.flat.W.YoY", "b.notflat.W.YoY", "n.b.notflat.W.YoY", "n.bU.W.YoY",
                 "SplineDesign.H.YoY",
                 "b.flat.H.YoY", "n.b.flat.H.YoY", "b.notflat.H.YoY", "n.b.notflat.H.YoY", "n.bU.H.YoY",
                 "SplineDesign.W.1",
                 "b.flat.W.1", "n.b.flat.W.1", "b.notflat.W.1", "n.b.notflat.W.1", "n.bU.W.1",
                 "SplineDesign.H.1",
                 "b.flat.H.1", "n.b.flat.H.1", "b.notflat.H.1", "n.b.notflat.H.1", "n.bU.H.1",
                 "tauU.alpha", "tauU.beta", "taueU.alpha", "taueU.beta",
                 "mu_xiP", "tau_xiP", "tauP.alpha", "tauP.beta")

parameters <- c("logitP", "beta.logitP", "tauP", "sigmaP",
                "bU.W.YoY", "bU.H.YoY", "tauU", "sigmaU", 
                "eU.W.YoY", "eU.H.YoY", "taueU", "sigmaeU", 
                "Utot.W.YoY", "Utot.H.YoY", "Utot.YoY", "logUne.W.YoY", "logUne.H.YoY",
                "etaU.W.YoY", "etaU.H.YoY", "U.W.YoY", "U.H.YoY",
                "bU.W.1", "bU.H.1",
                "eU.W.1", "eU.H.1", 
                "Utot.W.1", "Utot.H.1", "Utot.1", "logUne.W.1", "logUne.H.1",
                "etaU.W.1", "etaU.H.1", "U.W.1", "U.H.1",
                "Utot" )
if( any(is.na(m2))) {parameters <- c(parameters,"m2")} # monitor in case some bad data where missing values present
if( any(is.na(u2.A.YoY))) {parameters <- c(parameters,"u2.A.YoY")}
if( any(is.na(u2.N.YoY))) {parameters <- c(parameters,"u2.N.YoY")}
if( any(is.na(u2.A.1  ))) {parameters <- c(parameters,"u2.A.1")}
if( any(is.na(u2.N.1  ))) {parameters <- c(parameters,"u2.N.1")}


# Now to create the initial values, and the data prior to call to the MCMC sampler

Nstrata <- length(n1)

# Estimate number of YoY wild and hatchery fish based on clip rate
u2.H.YoY <- u2.A.YoY/clip.frac.H.YoY  # only a portion of the YoY hatchery fish are clipped
u2.W.YoY <- pmax(u2.N.YoY - u2.H.YoY*(1-clip.frac.H.YoY),0) # subtract the questimated number of hatchery fish
u2.H.YoY[is.na(u2.H.YoY)] <- 1  # in case of missing values
u2.W.YoY[is.na(u2.W.YoY)] <- 1  # in case of missing values
# Estimate number of Age1 wild and hatchery fish based on clip rate
u2.H.1   <- u2.A.1/clip.frac.H.1  # only a portion of the AGE1 hatchery fish are clipped
u2.W.1   <- pmax(u2.N.1 - u2.H.1*(1-clip.frac.H.1),0) # subtract the questimated number of hatchery fish
u2.H.1[is.na(u2.H.1)] <- 1  # in case of missing values
u2.W.1[is.na(u2.W.1)] <- 1  # in case of missing values

avg.P <- sum(m2,na.rm=TRUE)/sum(n1, na.rM=TRUE)
Uguess.W.YoY <- pmax((u2.W.YoY+1)*(n1+2)/(m2+1), u2.W.YoY/avg.P, 1, na.rm=TRUE)  # try and keep Uguess larger than observed values
Uguess.H.YoY <- pmax((u2.H.YoY+1)*(n1+2)/(m2+1), u2.H.YoY/avg.P, 1, na.rm=TRUE)
Uguess.H.YoY[1:hatch.after.YoY] <- 0   # no YoY hatchery fish prior to release from hatchery
Uguess.W.1   <- pmax((u2.W.1+1)*(n1+2)/(m2+1), u2.W.1/avg.P, 1, na.rm=TRUE)  # try and keep Uguess larger than observed values
Uguess.H.1   <- pmax((u2.H.1+1)*(n1+2)/(m2+1), u2.H.1/avg.P, 1, na.rm=TRUE)

# create the B-spline design matrix for YoY wild and hatchery fish
# The design matrix for hatchery fish will still have rows corresponding to entries PRIOR to 
#    the hatchery release but these are never used in the winbugs fitting routines
# There is a separate (single) spline for hatchery and wild fish with NO breakpoints
# The first two coefficient have a flat prior and the rest of the coefficients are modelled using
#    differences between the succesive coefficients

# YoY Wild fish. This covers the entire experiment.
SplineDegree <- 3           # Degree of spline between occasions
knots <- seq(4,Nstrata,4)/(Nstrata+1) # a knot roughly every 4th stratum
SplineDesign.W.YoY <- bs((1:Nstrata)/(Nstrata+1), knots=knots, degree=SplineDegree, intercept=TRUE, Boundary.knots=c(0,1))
b.flat.W.YoY      <- c(1,2)
b.notflat.W.YoY   <- 3:(ncol(SplineDesign.W.YoY))
n.b.flat.W.YoY    <- length(b.flat.W.YoY)
n.b.notflat.W.YoY <- length(b.notflat.W.YoY)
n.bU.W.YoY        <- n.b.flat.W.YoY + n.b.notflat.W.YoY
init.bU.W.YoY   <- lm(log(Uguess.W.YoY+1) ~ SplineDesign.W.YoY-1)$coefficients  # initial values for spline coefficients

# Age1 Wild fish. This covers the entire experiment.
SplineDegree <- 3           # Degree of spline between occasions
knots <- seq(4,Nstrata,4)/(Nstrata+1) # a knot roughly every 4th stratum
SplineDesign.W.1  <- bs((1:Nstrata)/(Nstrata+1), knots=knots, degree=SplineDegree, intercept=TRUE, Boundary.knots=c(0,1))
b.flat.W.1        <- c(1,2)
b.notflat.W.1     <- 3:(ncol(SplineDesign.W.1))
n.b.flat.W.1      <- length(b.flat.W.1)
n.b.notflat.W.1   <- length(b.notflat.W.1)
n.bU.W.1          <- n.b.flat.W.1 + n.b.notflat.W.1
init.bU.W.1     <- lm(log(Uguess.W.1+1) ~ SplineDesign.W.1-1)$coefficients  # initial values for spline coefficients

# YoY hatchery fish. Notice they can only enter AFTER hatch.after, The spline design matrix still has rows
# of zero for 1:hatch.after to make it easier in Bugs
SplineDegree <- 3           # Degree of spline between occasions
knots <- (seq((hatch.after.YoY+4),Nstrata-1,4)-hatch.after.YoY)/(Nstrata-hatch.after.YoY+1) # a knot roughly every 4th stratum
SplineDesign.H.YoY <- bs((1:(Nstrata-hatch.after.YoY))/(Nstrata-hatch.after.YoY+1), knots=knots, degree=SplineDegree, intercept=TRUE, Boundary.knots=c(0,1))
b.flat.H.YoY      <- c(1,2)
b.notflat.H.YoY   <- 3:(ncol(SplineDesign.H.YoY))
n.b.flat.H.YoY    <- length(b.flat.H.YoY)
n.b.notflat.H.YoY <- length(b.notflat.H.YoY)
n.bU.H.YoY        <- n.b.flat.H.YoY + n.b.notflat.H.YoY
init.bU.H.YoY   <- lm(log(Uguess.H.YoY[(hatch.after.YoY+1):Nstrata]+1) ~ SplineDesign.H.YoY-1)$coefficients  # initial values for spline coefficients
# patch up the initial rows of the spline design matrix
SplineDesign.H.YoY <- rbind(matrix(0,nrow=hatch.after.YoY, ncol=ncol(SplineDesign.H.YoY)), SplineDesign.H.YoY)

# Age1 hatchery fish. These are present from last year and so we must model the entire experiment
SplineDegree <- 3           # Degree of spline between occasions
knots <- seq(4,Nstrata,4)/(Nstrata+1) # a knot roughly every 4th stratum
SplineDesign.H.1 <- bs((1:Nstrata)/(Nstrata+1), knots=knots, degree=SplineDegree, intercept=TRUE, Boundary.knots=c(0,1))
b.flat.H.1        <- c(1,2)
b.notflat.H.1     <- 3:(ncol(SplineDesign.H.1))
n.b.flat.H.1      <- length(b.flat.H.1)
n.b.notflat.H.1   <- length(b.notflat.H.1)
n.bU.H.1          <- n.b.flat.H.1 + n.b.notflat.H.1
init.bU.H.1   <- lm(log(Uguess.H.1+1) ~ SplineDesign.H.1-1)$coefficients  # initial values for spline coefficients




# create an initial plot of the fit to the number of YoY and Age1 unmarked fish
pdf(file=paste(prefix,"-initialU.pdf",sep=""))
ylim <- c( min( c(log(Uguess.H.YoY+1),log(Uguess.W.YoY+1),log(Uguess.H.1+1),log(Uguess.W.1+1)), na.rm=TRUE),
           max( c(log(Uguess.H.YoY+1),log(Uguess.W.YoY+1),log(Uguess.H.1+1),log(Uguess.W.1+1)), na.rm=TRUE))
plot(time, log(Uguess.H.YoY+1), 
    main=paste(title,"\nInitial spline fit to estimated U.W[i] and U.H[i]"),
    sub="h,w = YoY, H,W=Age 1",
    ylab="log(U[i])", xlab='Stratum', pch="h", ylime=ylim)  # initial points on log scale.
points(time, log(Uguess.W.YoY+1), pch="w")
points(time, log(Uguess.H.1+1), pch="H")  # age1 fish
points(time, log(Uguess.W.1+1), pch="W")
lines(time, SplineDesign.W.YoY %*% init.bU.W.YoY)  # add smoothed spline through points
lines(time, SplineDesign.H.YoY %*% init.bU.H.YoY)  # add smoothed spline through points
lines(time, SplineDesign.W.1   %*% init.bU.W.1  )  # add smoothed spline through points
lines(time, SplineDesign.H.1   %*% init.bU.H.1  )  # add smoothed spline through points
dev.off()

#browser()

# get the logitP=logit(P) covariate matrix ready 
logitP.cov <- as.matrix(logitP.cov)
NlogitP.cov <- ncol(as.matrix(logitP.cov))


#  initial values for the parameters of the model
                
init.vals <- function(){
#  browser()
   # Initial values for the probability of capture
   init.logitP <- pmax(-10,pmin(10,logit((m2+1)/(n1+2))))         # initial capture rates based on observed recaptures
   init.logitP[is.na(init.logitP)] <- -2         # those cases where initial probability is unknown
   init.beta.logitP <- as.vector(lm( init.logitP ~ logitP.cov-1)$coefficients)
   init.beta.logitP[init.beta.logitP=NA] <- 0 
   init.beta.logitP <- c(init.beta.logitP, 0)   # add one extra element so that single beta is still written as a 
                                             # vector in the init files etc.
   init.tauP <- 1/var(init.logitP, na.rm=TRUE)     # 1/variance of logit(p)'s (ignoring the covariates for now)

   # inital values for the YoY spline coefficients
   init.bU.W.YoY   <- lm(log(Uguess.W.YoY+1) ~ SplineDesign.W.YoY-1)$coefficients  
   init.bU.H.YoY   <- lm(log(Uguess.H.YoY[(hatch.after.YoY+1):Nstrata]+1) ~ SplineDesign.H.YoY[(hatch.after.YoY+1):Nstrata,]-1)$coefficients  
   # inital values for the Age1 spline coefficients
   init.bU.W.1     <- lm(log(Uguess.W.1+1) ~ SplineDesign.W.1-1)$coefficients  
   init.bU.H.1     <- lm(log(Uguess.H.1+1) ~ SplineDesign.H.1-1)$coefficients  

   init.eU.W.YoY   <- as.vector(log(Uguess.W.YoY+1)-SplineDesign.W.YoY%*%init.bU.W.YoY)  # error terms set as differ between obs and pred
   init.etaU.W.YoY <- log(Uguess.W.YoY+1)
   init.eU.W.1     <- as.vector(log(Uguess.W.1  +1)-SplineDesign.W.1  %*%init.bU.W.1  )  # error terms set as differ between obs and pred
   init.etaU.W.1   <- log(Uguess.W.1  +1)

   init.eU.H.YoY   <- as.vector(log(Uguess.H.YoY+1)-SplineDesign.H.YoY%*%init.bU.H.YoY)  # error terms set as differ between obs and pred
   init.etaU.H.YoY <- log(Uguess.H.YoY+1)
#  init.etaU.H.YoY[1:hatch.after.YoY] <- NA  # these are never used.
   init.eU.H.1     <- as.vector(log(Uguess.H.1  +1)-SplineDesign.H.1  %*%init.bU.H.1  )  # error terms set as differ between obs and pred
   init.etaU.H.1   <- log(Uguess.H.1  +1)

   # variance of spline difference (use only the YoY wild fish to initialize)
   sigmaU <- sd( init.bU.W.YoY[b.notflat.W.YoY]-2*init.bU.W.YoY[b.notflat.W.YoY-1]+init.bU.W.YoY[b.notflat.W.YoY-2], na.rm=TRUE)
   init.tauU <- 1/sigmaU^2

   # variance of error in the U' over and above the spline fit (use only the YoY wild fish to initialize)
   sigmaeU <- sd(init.eU.W.YoY, na.rm=TRUE)
   init.taueU <- 1/sigmaeU^2

   # initialize the u2.A.YoY and u2.N.YoY where missing
   init.u2.A.YoY    <- u2.A.YoY
   init.u2.A.YoY[ is.na(u2.A.YoY)] <- 100
   init.u2.A.YoY[!is.na(u2.A.YoY)] <- NA
   init.u2.A.1      <- u2.A.1  
   init.u2.A.1  [ is.na(u2.A.1  )] <- 100
   init.u2.A.1  [!is.na(u2.A.1  )] <- NA

   init.u2.N.YoY    <- u2.N.YoY
   init.u2.N.YoY[ is.na(u2.N.YoY)] <- 100
   init.u2.N.YoY[!is.na(u2.N.YoY)] <- NA
   init.u2.N.1      <- u2.N.1  
   init.u2.N.1  [ is.na(u2.N.1  )] <- 100
   init.u2.N.1  [!is.na(u2.N.1  )] <- NA

   list(logitP=init.logitP, beta.logitP=init.beta.logitP, tauP=init.tauP, 
        bU.W.YoY=init.bU.W.YoY, bU.H.YoY=init.bU.H.YoY,   tauU=init.tauU, taueU=init.taueU, 
        etaU.W.YoY=init.etaU.W.YoY, etaU.H.YoY=init.etaU.H.YoY,
        bU.W.1  =init.bU.W.1  , bU.H.1  =init.bU.H.1  , 
        etaU.W.1=init.etaU.W.1, etaU.H.1=init.etaU.H.1)
}
# make a list of initial values
init.vals.list <- lapply(1:n.chains, function(x){init.vals()})

#browser()


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
                        initVals=init.vals.list,
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
