# 2014-09-01 CJS conversion to JAGS
#                - no model name
#                - C(,20) -> T(,20)
#                - dflat() to dnorm(0, 1E-6)
#                - added u2copy to improve mixing based on Matt S. suggestion
# 2011-05-15 CJS Limited etaU to max of 20
# 2011-02-28 CJS First version

TimeStratPetersenNonDiagErrorNPMarkAvail <- function(title,
                                                     prefix,
                                                     time,
                                                     n1,
                                                     m2,
                                                     u2,
                                                     jump.after=NULL,
                                                     logitP.cov=rep(1,length(u2)),
                                                     logitP.fixed=rep(NA,length(u2)),
                                                     ma.p.alpha,
                                                     ma.p.beta,
                                                     n.chains=3,
                                                     n.iter=200000,
                                                     n.burnin=100000,
                                                     n.sims=2000, 
                                                     tauU.alpha=1,
                                                     tauU.beta=.05,
                                                     taueU.alpha=1,
                                                     taueU.beta=.05,
                                                     Delta.max,
                                                     tauTT.alpha=.1,
                                                     tauTT.beta=.1,
                                                     mu_xiP=-2,
                                                     tau_xiP=.6666, 
                                                     tauP.alpha=.001,
                                                     tauP.beta=.001, 
                                                     debug=FALSE,
                                                     debug2=FALSE,
						     engine=c('jags','openbugs')[1],
                                                     InitialSeed){

set.seed(InitialSeed)  # set prior to initial value computations

#
#  Fit the smoothed time-Stratified Petersen estimator with NON-Diagonal recoveries and less than 100%
#  marked fish available for recapture.
#  This model allows recoveries outside the stratum of release and error in the smoothed U curve.
#  The travel time model is based on the continuation ratio and makes no parametric assumptions.
#  It allows for only a fraction of marked fish to be available based on prior information
#  on the marked availability rate (ma.p).
#
#  This routine assumes that the strata are time (e.g. weeks).
#  In each stratum n1 fish are released (with marks). These are ususally
#     captured fish that are marked, transported upstream, and released.
#     These fish are used only to estimate the recapture rate downstream.
#  Not all of the marked fish are available for subsequent recapture. Only a fraction ma.p
#     are available. This lack of availability could be because of fall back, handling mortality, etc.
#  Of the n1 fish released and available, some are recapturd in this stratum of release (column 1 of m2) or in
#     subsequent strata (subsequent columns of m2). No fish are assumed to be available for capture
#     outside the range of strata considered in the matrix of m2
#  At the same tine, u2 other (unmarked) fish are newly captured in stratum i. 
#     These EXCLUDE recaptures of marked fish. These are the fish that are "expanded"
#     to estimate the population size of fish in stratum i.  
# 
#  Input
#      prefix - prefix for file name for initial plot of U's
#      time- the stratum number
#      n1  - vector of number of fish released in stratum i
#      m2  - matrix of number of fish recovered who were released in stratum i and recovered in stratum j
#      u2  - vector of number of unmarked fish captured in stratum i
#      jump.after - points after which the spline is allowed to jump. Specify as a list of integers in the 
#              range of 1:Nstrata. If jump.after[i]=k, then the spline is split between strata k and k+1
#      logitP.cov - covariates for logit(P)=X beta.logitP.cov
#                 - specify anything you want for fixed logitP's as the covariate values are simply ignored. 
#                 - recommend that you specify 1 for the intercept and 0's for everything else
#      logitP.fixed - values for logitP that are fixed in advance. Use NA if corresponding value is not fixed, 
#                    otherwise specify the logitP value.
#      ma.p.alpha, ma.p.beta - information on mark available. Assumed to be prior beta(ma.p.alpha, ma.p.beta)


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

# Time Stratified Petersen with NON Diagonal recapture and allowing for error in the smoothed U curve.
# Non-parametric estimateion of travel times for marked individuals.
#
# Refer to Bonner and Schwarz (2010) Smoothed Estimates
#   for Time-Stratified Mark-Recapture Experiments using Bayesian
#   P-Splines
#
#  Data input:
#      Nstrata.rel - number of strata where fish are releases
#      Nstrata.cap - number of (future strata) where fish are recaptured.
#      n1         - number of marked fish released
#      m2         - number of marked fish recaptured
#                   This is a matrix of size Nstrata.rel x (Nstrata.cap+1)
#                   with entries m2[i,j] = number of fish released in i and recaptured in j
#                   Entries in the last column are the number of fish NEVER recaptured from those
#                   released
#      u2         - number of unmarked fish captured (To be expanded to population).
#      logitP     - the recapture rates. Use NA if these are modelled, otherwise specify the logit(fixed value, e.g. -10 for 0).
#      Nfree.logitP - number of free logitP parameters
#      free.logitP.index - vector of length(Nfree.logitP) for the free logitP parameters
#      logitP.cov - covariates for logitP
#      NlogitP.cov  - number of logitP covariates
#      ma.p.alpha, ma.p.beta - prior beta(alpha,beta) on mark availability
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
#                       - the other beta terms are given a prior of a N(mu=0, variance=1000)
#      tauP.alpha, tauP.beta - parameter for prior on tauP (residual variance of logit(P)'s after adjusting for
#                         covariates)
#      xiMu, tauMu  - mean and precision (1/variance) for prior on mean(log travel-times)
#      siSd, tauSd  - mean and precision (1/variance) for prior on sd(log travel times)
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
#
#      muTT[j] = mean(logit(delta[i,i+j-1])), j=1,...,Delta.max
#      sdTT = sd(logit(delta[i,i+j-1])), j=1,....,Delta.max
#      delta[i,i+j-1]=Theta[i,i+j-1]/(1-Theta[i,i]-...-Theta[i,i+j-2])
#       

   ##### Fit the spline for the U's and specify hierarchial model for the logit(P)'s ######
   for(i in 1:(Nstrata.cap)){
        ## Model for U's
        logUne[i] <- inprod(SplineDesign[i,1:n.bU],bU[1:n.bU])  # spline design matrix * spline coeff 
", fill=TRUE)
sink()  # Temporary end of saving bugs program
if(tolower(engine)=="jags") {
   sink("model.txt", append=TRUE)
   cat("
        etaU[i] ~ dnorm(logUne[i], taueU)T(,20)              # add random error
   ",fill=TRUE)
   sink()
}
if(tolower(engine) %in% c("openbugs")) {
   sink("model.txt", append=TRUE)
   cat("
        etaU[i] ~ dnorm(logUne[i], taueU)C(,20)              # add random error but keep from getting too large
   ",fill=TRUE)
   sink()
}
   sink("model.txt", append=TRUE)
   cat("
        eU[i] <- etaU[i] - logUne[i]
   }

   for(i in 1:Nfree.logitP){   # model the free capture rates using covariates
        mu.logitP[free.logitP.index[i]] <- inprod(logitP.cov[free.logitP.index[i],1:NlogitP.cov], beta.logitP[1:NlogitP.cov]) 
        logitP[free.logitP.index[i]] ~ dnorm(mu.logitP[free.logitP.index[i]],tauP)
   }

   ##### Priors and hyperpriors #####

   ##### Prior information on mark availability #####
   ##    There is no other information in the actual study to update the ma.p other than the prior
   ##    information.
   ma.p ~ dbeta(ma.p.alpha, ma.p.beta)

   ## Transition probabilities -- continuation ratio model
   for(i in 1:Nstrata.rel){
        ## delta[i,j] is the probability that a marked fish released on day i passes the second trap
        ## on day i+j-1 given that it does not pass the on days i,...,i+j-2. r[i,j]=logit(delta[i,j])
        ## is assumed to have a normal distribution with mean muTT[j] and precision tauTT.

        r[i,1] ~ dnorm(muTT[1],tauTT)
		
	logit(Theta[i,1] ) <- r[i,1]
		
	for(j in 2:Delta.max){
		r[i,j] ~ dnorm(muTT[j],tauTT)
			
		logit(delta[i,j]) <- r[i,j]
			
		Theta[i,j] <- delta[i,j] * (1 - sum(Theta[i,1:(j-1)]))
	}
	Theta[i,Delta.max+1] <- 1- sum(Theta[i,1:Delta.max])
   }

   for(j in 1:Delta.max){
	muTT[j] ~ dnorm(0,.666)
   }
   tauTT~ dgamma(tauTT.alpha,tauTT.beta)
   sdTT <- 1/sqrt(tauTT)
	
   ## Run size - flat priors
      ", fill=TRUE)
sink()  # Temporary end of saving bugs program
if(tolower(engine)=="jags") {
   sink("model.txt", append=TRUE)
   cat("
   for(i in 1:n.b.flat){
      bU[b.flat[i]] ~ dnorm(0, 1E-6)
   }
   ",fill=TRUE)
   sink()
}
if(tolower(engine) %in% c("openbugs")) {
   sink("model.txt", append=TRUE)
   cat("
   for(i in 1:n.b.flat){
      bU[b.flat[i]] ~ dflat()
   }
   ",fill=TRUE)
   sink()
}
   sink("model.txt", append=TRUE)
   cat("

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
   ## marked fish ##
   for(i in 1:Nstrata.rel){
      # Compute cell probabilities
      for(j in 1:(Delta.max+1)){
        Pmarked[i,j] <- Theta[i,j] * p[i+j-1] * ma.p  # adjust for availability
      }
      Pmarked[i,Delta.max+2] <- 1- sum(Pmarked[i,1:(Delta.max+1)])

      # Likelihood contribution
      m2[i,1:(Delta.max+2)] ~ dmulti(Pmarked[i,],n1[i])
   }

   ## Capture probabilities and run size
   for(j in 1:(Nstrata.cap + Extra.strata.cap)){
      logit(p[j]) <- logitP[j]       # convert from logit scale
   }
   for(j in 1:Nstrata.cap){
      U[j]   <- round(exp(etaU[j]))       # convert from log scale
      u2[j] ~ dbin(p[j],U[j])      # capture of newly unmarked fish
   }

   ##### Derived Parameters #####
   Utot <- sum( U[1:Nstrata.cap])          # Total number of unmarked fish
   Ntot <- sum(n1[1:Nstrata.rel]) + Utot  # Total population size including those fish marked and released
} # end of model
", fill=TRUE)

sink()  # End of saving the Bugs program

# Now to create the initial values, and the data prior to call to the MCMC sampler

Nstrata.rel <- length(n1)
Nstrata.cap <- length(u2)

## Count extra columns that will have to be added to account for Delta.max
Extra.strata.cap <- max(0,Nstrata.rel + ncol(m2) - Nstrata.cap -1)

# create the B-spline design matrix
# Each set of strata separated at the jump.after[i] points forms a separate spline with a separate basis
# We need to keep track of the breaks as the first two spline coefficients will have a flat
# prior and the others are then related to the previous values.

ext.jump <- c(0, jump.after, Nstrata.cap)  # add the first and last breakpoints to the jump sets
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


# get the logitP covariate matrix ready 
logitP.cov <- as.matrix(logitP.cov)
NlogitP.cov <- ncol(as.matrix(logitP.cov))

# get the logitP's ready to allow for fixed values
logitP <- c(as.numeric(logitP.fixed),rep(-10,Extra.strata.cap))
storage.mode(logitP) <- "double"  # force the storage class to be correct if there are no fixed values
free.logitP.index <- (1:Nstrata.cap)[ is.na(logitP.fixed)]  # free values are those where NA is specifed
Nfree.logitP <- length(free.logitP.index)

# make a copy of u2 to improve mixing (not yet implemented)
#u2copy <- spline(x=1:length(u2), y=u2, xout=1:length(u2))$y
#u2copy <- round(u2copy) # round to integers

datalist <- list("Nstrata.rel", "Nstrata.cap","Extra.strata.cap",
                 "Delta.max","n1", "m2", "u2", # "u2copy", # u2copy not yet implemented
                 "logitP", "Nfree.logitP", "free.logitP.index",
                 "logitP.cov", "NlogitP.cov",
                 "ma.p.alpha","ma.p.beta",
                 "SplineDesign",
                 "b.flat", "n.b.flat", "b.notflat", "n.b.notflat", "n.bU",
                 "tauTT.alpha","tauTT.beta",
                 "tauU.alpha", "tauU.beta", "taueU.alpha", "taueU.beta",
                 "mu_xiP", "tau_xiP", "tauP.alpha", "tauP.beta")


## Generate the initial values for the parameters of the model

## 1) U and spline coefficients
Uguess <- pmax((u2+1)/expit(mu_xiP),1)  # try and keep Uguess larger than observed values
Uguess[which(is.na(Uguess))] <- mean(Uguess,na.rm=TRUE)

init.bU   <- lm(log(Uguess) ~ SplineDesign-1)$coefficients  # initial values for spline coefficients

if(debug2) {
   cat("compute init.bU \n")
   browser()  # Stop here to examine the spline design matrix function
}

## 2) Capture probabilities
logitPguess <- c(logit(pmin(.99,pmax(.01,(apply(m2[,1:(Delta.max+1)],1,sum)+1)/(n1+1)))),
                 rep(mu_xiP,Nstrata.cap-Nstrata.rel))
#browser()
init.beta.logitP <- as.vector(lm( logitPguess ~ logitP.cov-1)$coefficients)
if(debug2) {
   cat(" obtained initial values of beta.logitP\n")
   browser()
}

## 3) marked availability
ma.p.guess <- ma.p.alpha/(ma.p.alpha+ma.p.beta)


# create an initial plot of the fit
pdf(file=paste(prefix,"-initialU.pdf",sep=""))
plot(time, log(Uguess[1:Nstrata.cap]), 
    main=paste(title,"\nInitial spline fit to estimated U[i]"),
    ylab="log(U[i])", xlab='Stratum')  # initial points on log scale.
lines(time, SplineDesign %*% init.bU)  # add smoothed spline through points
dev.off()

parameters <- c("logitP", "beta.logitP", "tauP", "sigmaP",
                "bU", "tauU", "sigmaU", 
                "eU", "taueU", "sigmaeU", 
                "Ntot", "Utot", "logUne", "etaU", "U",
                 "muTT","sdTT","Theta","ma.p")

if( any(is.na(m2))) {parameters <- c(parameters,"m2")} # monitor in case some bad data where missing values present
if( any(is.na(u2))) {parameters <- c(parameters,"u2")}
                 
init.vals <- function(){
   init.logitP <- c(logit((apply(m2[,1:(Delta.max+1)],1,sum)+1)/(n1+1)),rep(mu_xiP,Nstrata.cap-Nstrata.rel))         # initial capture rates based on observed recaptures
   init.logitP <- pmin(10,pmax(-10,init.logitP))
   init.logitP[is.na(init.logitP)] <- -2         # those cases where initial probability is unknown
   init.logitP[!is.na(logitP.fixed)] <- NA        # no need to initialize the fixed values
   
   init.beta.logitP <- as.vector(lm( init.logitP ~ logitP.cov-1)$coefficients)
   init.beta.logitP[is.na(init.beta.logitP)] <- 0 
   init.beta.logitP <- c(init.beta.logitP, 0)   # add one extra element so that single beta is still written as a 
                                             # vector in the init files etc.

   init.logitP <- c(init.logitP,rep(NA,Extra.strata.cap)) # Add values for extra capture probabilities
   
   init.tauP <- 1/var(init.logitP, na.rm=TRUE)     # 1/variance of logit(p)'s (ignoring the covariates for now)

   init.bU   <- lm(log(Uguess) ~ SplineDesign-1)$coefficients  # initial values for spline coefficients
   init.eU   <- as.vector(log(Uguess)-SplineDesign%*%init.bU)  # error terms set as differ between obs and pred
   init.etaU <- log(Uguess)

   # variance of spline difference
   sigmaU <- sd( init.bU[b.notflat]-2*init.bU[b.notflat-1]+init.bU[b.notflat-2], na.rm=TRUE)
   init.tauU <- 1/sigmaU^2

   # variance of error in the U' over and above the spline fit
   sigmaeU <- sd(init.eU, na.rm=TRUE)
   init.taueU <- 1/sigmaeU^2

   # initialize the u2 where missing
   init.u2    <- u2
   init.u2[ is.na(u2)] <- 100
   init.u2[!is.na(u2)] <- NA

   ## Transition probabilities
   init.Theta <- t(sapply(1:Nstrata.rel,function(i){

     if(all(is.na(m2[i,])) || sum(m2[i,])==0)
        return(rep(NA,Delta.max+1))

     else{
       thetatmp <- pmax(.01,pmin(m2[i,-(Delta.max+2)]/sum(m2[i,-(Delta.max+2)],na.rm=TRUE),.99,na.rm=TRUE))  # CJS 2011-02-16
       return(thetatmp/sum(thetatmp))
     }
   }))
# cat('Initial values')
# browser()
   init.delta <- t(apply(as.matrix(init.Theta[,-(Delta.max+1)]),1,   # CJS 2011-02-16 as.matrix added
      function(theta){    # CJS fixed -(Delta.max+1)
         if(length(theta) == 1){theta} else {theta/(1-c(0,cumsum(theta[-Delta.max])))}
   }))
   
   ## mean and standard deviation of transition probabilties
   init.muTT <- apply(logit(init.delta),2,mean,na.rm=TRUE)
   init.sdTT <- sd(as.vector(t(logit(init.delta)))-init.muTT,na.rm=TRUE)
 
   ## ma.p
   init.ma.p <- ma.p.alpha/(ma.p.alpha+ma.p.beta)
  
#browser() 
   list(logitP=init.logitP, beta.logitP=init.beta.logitP, tauP=init.tauP, 
        bU=init.bU,  tauU=init.tauU, taueU=init.taueU, etaU=init.etaU, 
        muTT=init.muTT, tauTT=1/init.sdTT^2,r=logit(init.delta), ma.p=init.ma.p)
}

#browser()


## Generate data list
data.list <- list()
for(i in 1:length(datalist)){
  data.list[[length(data.list)+1]] <-get(datalist[[i]])
}
names(data.list) <- as.list(datalist)

## Generate the initial values and put into a list
# make a list of initial values
init.vals.list <- lapply(1:n.chains, function(x){init.vals()})
  
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
