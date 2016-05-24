# 2014-09-01 CJS Converted to JAGS
## 2013-12-31 CJS Tried adding u2copy to get back Matts fix for mixing
## 2013-09-22 SJB Changes to model for JAGS compatability:
##     -- Removed model name.
##     -- Changed C(,20) to T(,20).
##     -- Replace dflat() with dnorm(0.0,1.0E-6).
##     -- Remove Matt's fix to improve mixing.
# 2013-09-04 CJS removed all references to WinBugs. Fixed problem with initial values for NA in n1, m2, or u2
# 2011-05-15 CJS limited etaU to 20 or less
# 2011-01-24 SB  added call to run.windows.openbugs and run.windows.winbugs
# 2010-11-25 CJS add output to show progress of sampling through burnin and post-burnin phases
# 2010-04-26 CJS fixed problem in computing logitPguess when m2=n1 and you get infinite logit value
# 2009-12-05 CJS added title to argument list
# 2009-12-01 CJS (added WinBugs/OpenBugs directory to the argument list

TimeStratPetersenDiagError <- function(
    title,
    prefix,
    time,
    n1,
    m2,
    u2,
    jump.after=NULL,
    logitP.cov,
    n.chains=3,
    n.iter=200000,
    n.burnin=100000,
    n.sims=2000,
    tauU.alpha=1,
    tauU.beta=.05,
    taueU.alpha=1,
    taueU.beta=.05,
    mu_xiP=logit(sum(m2,na.rm=TRUE)/sum(n1,na.rm=TRUE)),
    tau_xiP=1/var(logit((m2+.5)/(n1+1)),na.rm=TRUE),
    tauP.alpha=.001, tauP.beta=.001,
    debug=FALSE,
    debug2=FALSE,
    engine=c("jags","openbugs")[1],
    InitialSeed){

set.seed(InitialSeed)  # set prior to initial value computations

#
#  Fit the smoothed time-Stratified Petersen estimator with Diagonal recoveries (i.e. no recoveries
#  outside stratum of release) and error in the smoothed U curve
#
#  Packages Required - must be installed BEFORE calling this functin
#
#    R2OpenBugs  - needed for the as.bugs.array() function
#    BRugs       - only needed if using OpenBugs
#    rjags       - only needed if using JAGS
#    Coda
#    actuar
#    splines
#
#  This routine assumes that the strata are time (e.g. weeks).
#  In each stratum n1 fish are released (with marks). These are ususall
#     captured fish that are marked, transported upstream, and released.
#     These fish are used only to estimate the recapture rate downstream.
#  Of the n1 fish released, m2 fish are recaptured in the same stratum (e.g. week) of release.
#     There is a related function that allows fish to be recaptured in subsequent weeks.
#  At the same tine, u2 other (unmarked) fish are newly captured in stratum i.
#     These EXCLUDE recaptures of marked fish. These are the fish that are "expanded"
#     to estimate the population size of fish in stratum i.
#
#  Input
#      prefix - prefix for file name for initial plot of U's
#      time- the stratum number
#      n1  - vector of number of fish released in stratum i
#      m2  - vector of number of fish recovered in stratum i (EXCLUDING recaps)
#      u2  - vector of number of unmarked fish captured in stratum i
#      jump.after - points after which the spline is allowed to jump. Specify as a list of integers in the
#              range of 1:Nstrata. If jump.after[i]=k, then the spline is split between strata k and k+1
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
model{
# Time Stratified Petersen with Diagonal recapture (no spillover in subsequent weeks or marked fish)
#    and allowing for error in the smoothed U curve.

# Refer to Bonner (2008) Ph.D. thesis from Simon Fraser University available at
#     http://www.stat.sfu.ca/people/alumni/Theses/Bonner-2008.pdf
# The model is in Appendix B. The discussion of the model is in Chapter 2.

#  Data input:
#      Nstrata - number of strata
#      n1         - number of marked fish released
#      m2         - number of marked fish recaptured
#      u2         - number of unmarked fish captured (To be expanded to population).
#      logitP.cov   - covariates for logitP
#      NlogitP.cov  - number of logitP covariates
#      SplineDesign- spline design matrix of size [Nstrata, maxelement of n.b.notflat]
#                   This is set up prior to the call.
#      b.flat   - vector of strata indices where the prior for the b's will be flat.
#                 this is normally the first two of each spline segment
#      n.b.flat - number of b coefficients that have a flat prior
#      b.notflat- vector of strata indices where difference in coefficients is modelled
#      n.b.notflat- number of b coefficients that do not have a flat prior
#      tauU.alpha, tauU.beta - parameters for prior on tauU
#      taueU.alpha, taueU.beta - parameters for prior on taueU
#      mu_xiP, tau_xiP  - parameters for prior on mean logit(P)'s [The intercept term]
#                         mu_xiP = mean; tau_xiP = 1/variance
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
#      U[i]
#       etaU[i]  = log(U[i])
#         which comes from spline with parameters bU[1... Knots+q]
#         + error term eU[i]

   ##### Fit the spline and specify hierarchial model for the logit(P)'s ######
   for(i in 1:Nstrata){
        logUne[i] <- inprod(SplineDesign[i,1:n.bU],bU[1:n.bU])  # spline design matrix * spline coeff
", fill=TRUE)
sink()  # Temporary end of saving bugs program
if(tolower(engine)=="jags") {
   sink("model.txt", append=TRUE)
   cat("
        etaU[i] ~ dnorm(logUne[i], taueU)T(,20)    # add random error
   ",fill=TRUE)
   sink()
}
if(tolower(engine) %in% c("openbugs")) {
   sink("model.txt", append=TRUE)
   cat("
        etaU[i] ~ dnorm(logUne[i], taueU)C(,20)    # add random error
   ",fill=TRUE)
   sink()
}
   sink("model.txt", append=TRUE)
   cat("
        eU[i] <- etaU[i] - logUne[i]
        mu.logitP[i] <- inprod(logitP.cov[i,1:NlogitP.cov], beta.logitP[1:NlogitP.cov])

#        logitPu[i] ~ dnorm(mu.logitP[i],tauP)       # uncontrained logitP value
#        logitP [i] <- max(-10,min(10, logitPu[i]))  # keep the logit from wandering too far off

        ## Matt's fix to improve mixing. Use u2copy to break the cycle (this doesn't work??)
        mu.epsilon[i] <- mu.logitP[i] - log(u2copy[i] + 1) + etaU[i]   
        epsilon[i] ~ dnorm(mu.epsilon[i],tauP)                     
        logitP[i] <- max(-10,min(10,log(u2copy[i] + 1) - etaU[i] + epsilon[i]))  
   }

   ##### Hyperpriors #####
   ## Run size - flat priors
   for(i in 1:n.b.flat){
", fill=TRUE)
sink()  # Temporary end of saving bugs program
if(tolower(engine)=="jags") {
   sink("model.txt", append=TRUE)
   cat("
      bU[b.flat[i]] ~ dnorm(0.0,1.0E-6) 
   ",fill=TRUE)
   sink()
}
if(tolower(engine) %in% c("openbugs")) {
   sink("model.txt", append=TRUE)
   cat("
      bU[b.flat[i]] ~ dflat()
   ",fill=TRUE)
   sink()
}
   sink("model.txt", append=TRUE)
   cat("
   }
   ## Run size - priors on the difference
   for(i in 1:n.b.notflat){
      xiU[b.notflat[i]] <- 2*bU[b.notflat[i]-1] - bU[b.notflat[i]-2]
      bU [b.notflat[i]] ~ dnorm(xiU[b.notflat[i]],tauU)
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
   for(i in 1:Nstrata){
      logit(p[i]) <- logitP[i]       # convert from logit scale
      U[i]   <- round(exp(etaU[i]))       # convert from log scale
      m2[i] ~ dbin(p[i],n1[i])     # recovery of marked fish
      u2[i] ~ dbin(p[i],U [i])      # capture of newly unmarked fish
   }

   ##### Derived Parameters #####
   Utot <- sum( U[1:Nstrata])          # Total number of unmarked fish
   Ntot <- sum(n1[1:Nstrata]) + Utot  # Total population size including those fish marked and released
} # end of model

", fill=TRUE)
sink()  # End of saving the Bugs program


# create the B-spline design matrix
# Each set of strata separated at the jump.after[i] points forms a separate spline with a separate basis
# We need to keep track of the breaks as the first two spline coefficients will have a flat
# prior and the others are then related to the previous values.

Nstrata <- length(n1)
ext.jump <- c(0, jump.after, Nstrata)  # add the first and last breakpoints to the jump sets
SplineDesign <- matrix(0, nrow=0, ncol=0)
SplineDegree <- 3           # Degree of spline between occasions
b.flat <- NULL              # index of spline coefficients with a flat prior distribution -first two of each segment
b.notflat <- NULL           # index of spline coefficients where difference is modelled
all.knots <- NULL
for (i in 1:(length(ext.jump)-1)){
  nstrata.in.set <- ext.jump[i+1]-ext.jump[i]
  if(nstrata.in.set > 7)
    { knots   <- seq(5,nstrata.in.set-1,4)/(nstrata.in.set+1) # a knot roughly every 4th stratum
    } else{
      knots   <- .5       # a knot roughly every 4th stratum
    }
  all.knots <- c(all.knots, knots)
  # compute the design matrix for this set of strata
  z <- bs((1:nstrata.in.set)/(nstrata.in.set+1), knots=knots, degree=SplineDegree,
             intercept=TRUE, Boundary.knots=c(0,1))
  # first two elements of b coeffients have a flat prior
  b.flat <- c(b.flat, ncol(SplineDesign)+(1:2))
  b.notflat <- c(b.notflat, ncol(SplineDesign)+3:(ncol(z)))
  # add to the full design matrix which is block diagonal
  SplineDesign <- cbind(SplineDesign, matrix(0, nrow=nrow(SplineDesign), ncol=ncol(z)))
  SplineDesign <- rbind(SplineDesign,
                         cbind( matrix(0,nrow=nrow(z),ncol=ncol(SplineDesign)-ncol(z)), z)  )
  } # end of for loop
n.b.flat <- length(b.flat)
n.b.notflat <- length(b.notflat)
n.bU <- n.b.flat + n.b.notflat


# get the logitP=logit(P) covariate matrix ready
logitP.cov <- as.matrix(logitP.cov)
NlogitP.cov <- ncol(as.matrix(logitP.cov))

# create a copy of the u2 to improve mixing in the MCMC model
u2copy <- spline(x=1:Nstrata, y=u2, xout=1:Nstrata)$y
u2copy <- round(u2copy)  # round to integers

datalist <- list("Nstrata", "n1", "m2", "u2", "u2copy", 
		 "logitP.cov", "NlogitP.cov",
                 "SplineDesign",
                 "b.flat", "n.b.flat", "b.notflat", "n.b.notflat", "n.bU",
                 "tauU.alpha", "tauU.beta", "taueU.alpha", "taueU.beta",
                 "mu_xiP", "tau_xiP", "tauP.alpha", "tauP.beta")


## Generate best guess initial values
## These initial values are used only to draw an initial fitted plot
## and are not used as initial values in the MCMC.

avgP <- sum(m2,na.rm=TRUE)/sum(n1,na.rm=TRUE)
Uguess <- pmax((u2+1)*(n1+2)/(m2+1), u2/avgP, 1, na.rm=TRUE)  # try and keep Uguess larger than observed values
Uguess[which(is.na(Uguess))] <- mean(Uguess,na.rm=TRUE)

init.bU   <- lm(log(Uguess+1) ~ SplineDesign-1)$coefficients  # initial values for spline coefficients
if(debug2) {
   cat("compute init.bU \n")
   browser()  # Stop here to examine the spline design matrix function
}

logitPguess <- pmax(-10, pmin(10, logit( (m2+1)/(n1+1))))
init.beta.logitP <- as.vector(lm( logitPguess ~ logitP.cov-1)$coefficients)
if(debug2) {
   cat(" obtained initial values of beta.logitP\n")
   browser()
}


# create an initial plot of the fit
pdf(file=paste(prefix,"-initialU.pdf",sep=""))
plot(time, log(Uguess),
    main=paste(title,"\nInitial spline fit to estimated U[i]"),
    ylab="log(U[i])", xlab='Stratum')  # initial points on log scale.
lines(time, SplineDesign %*% init.bU)  # add smoothed spline through points
dev.off()


parameters <- c("logitP", "beta.logitP", "tauP", "sigmaP",
                "bU", "tauU", "sigmaU",
                "eU", "taueU", "sigmaeU",
                "Ntot", "Utot", "logUne", "etaU", "U")
if( any(is.na(m2))) {parameters <- c(parameters,"m2")} # monitor in case some bad data where missing values present
if( any(is.na(u2))) {parameters <- c(parameters,"u2")}

## init.vals <- function(){
##    init.logitP <- logit((m2+1)/(n1+2))         # initial capture rates based on observed recaptures
##    init.logitP[is.na(init.logitP)] <- -2         # those cases where initial probability is unknown
##    init.beta.logitP <- as.vector(lm( init.logitP ~ logitP.cov-1)$coefficients)
##    init.beta.logitP[is.na(init.beta.logitP)] <- 0
##    init.beta.logitP <- c(init.beta.logitP, 0)   # add one extra element so that single beta is still written as a
##                                              # vector in the init files etc.
##    init.tauP <- 1/var(init.logitP, na.rm=TRUE)     # 1/variance of logit(p)'s (ignoring the covariates for now)

##    init.bU   <- lm(log(Uguess+1) ~ SplineDesign-1)$coefficients  # initial values for spline coefficients
##    init.eU   <- as.vector(log(Uguess)-SplineDesign%*%init.bU)  # error terms set as differ between obs and pred
##    init.etaU <- log(Uguess)

##    # variance of spline difference
##    sigmaU <- sd( init.bU[b.notflat]-2*init.bU[b.notflat-1]+init.bU[b.notflat-2], na.rm=TRUE)
##    init.tauU <- 1/sigmaU^2

##    # variance of error in the U' over and above the spline fit
##    sigmaeU <- sd(init.eU, na.rm=TRUE)
##    init.taueU <- 1/sigmaeU^2

##    # initialize the u2 where missing
##    init.u2    <- u2
##    init.u2[ is.na(u2)] <- 100
##    init.u2[!is.na(u2)] <- NA

##    list(logitP=init.logitP, beta.logitP=init.beta.logitP, tauP=init.tauP,
##         bU=init.bU,  tauU=init.tauU, taueU=init.taueU, etaU=init.etaU)
## }

## Generate initial values
init.vals <- genInitVals(model="TSPDE",
                         n1=n1,
                         m2=m2,
                         u2=u2,
                         logitP.cov=logitP.cov,
                         SplineDesign=SplineDesign,
                         n.chains=n.chains)

## Generate data list
data.list <- list()
for(i in 1:length(datalist)){
  data.list[[length(data.list)+1]] <-get(datalist[[i]])
}
names(data.list) <- as.list(datalist)

# Make the call to the MCMC sampler

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
