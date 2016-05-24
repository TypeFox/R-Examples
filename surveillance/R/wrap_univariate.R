##############################################################################
# This function is a wrapper for univariate surveillance algorithms
# using the old disProg and survRes object
#
# An sts object is given and a pre specified algorithms is ran
# by successively creating a disProg object for each region,
# running the algo and then assign the slots of the resulting survRes
# object to an sts object.
###################################################################################

###Apply other algorithms by wrapping up a suitable package.

#Wrapper function to call algo.farrington for each time series in an sts object
wrap.algo <- function(sts, algo, control,
                      control.hook=function(k, control) return(control),
                      verbose=TRUE,...) {
  #Number of time series
  nAreas <- ncol(sts@observed)
  nTimePoints <- nrow(sts@observed)
  nAlarm <- length(control$range)

  #Create alarm matrix having same size as sts
  sts@alarm <- matrix(NA,ncol=nAreas,nrow=nTimePoints,dimnames=dimnames(sts@observed))
  sts@upperbound <- matrix(NA,ncol=nAreas,nrow=nTimePoints,dimnames=dimnames(sts@observed))

  #Loop over all regions
  for (k in 1:nAreas) {
    if (verbose) {
      cat("Running ",algo," on area ",k," out of ",nAreas,"\n")
    }

    ##Create an old S3 disProg object
    disProg.k <- sts2disProg(sts[,k])

    #Use the univariate algorithm (possibly preprocess control object)
    kcontrol <- control.hook(k, control)
    survRes.k <- do.call(algo,args = list(disProg.k, control=kcontrol))

    #Transfer results to the S4 object
    if (!is.null(survRes.k)) {
      sts@alarm[control$range,k] <- survRes.k$alarm
      sts@upperbound[control$range,k] <- survRes.k$upperbound

      #Control object needs only to be set once
      sts@control <- survRes.k$control
    }
  }

  #Reduce sts object to only those obervations in range
  sts@observed <- sts@observed[control$range,,drop=FALSE]
  sts@state <- sts@state[control$range,,drop=FALSE]
  sts@populationFrac <- sts@populationFrac[control$range,,drop=FALSE]
  sts@alarm <- sts@alarm[control$range,,drop=FALSE]
  sts@upperbound <- sts@upperbound[control$range,,drop=FALSE]

  #Set correct theta0t matrix for all
  sts@control$theta0t <- control$theta0t

  #Fix the corresponding start entry
  start <- sts@start
  new.sampleNo <- start[2] + min(control$range) - 1
  start.year <- start[1] + (new.sampleNo - 1) %/% sts@freq
  start.sampleNo <- (new.sampleNo - 1) %% sts@freq + 1
  sts@start <- c(start.year,start.sampleNo)
  sts@epoch <- sts@epoch[control$range]
  sts@epochAsDate <- sts@epochAsDate

  #Ensure dimnames in the new object
  sts <- fix.dimnames(sts)

  return(sts)
}

#Farrington wrapper
farrington <- function(sts, control=list(range=NULL, b=3, w=3, reweight=TRUE, verbose=FALSE,alpha=0.01),...) {
  wrap.algo(sts,algo="algo.farrington",control=control,...)
}

#Bayes wrapper (this can be implemented more efficiently)
bayes <- function(sts, control = list(range = range, b = 0, w = 6, actY = TRUE,alpha=0.05),...) {
  if (sts@epochAsDate) {
    warning("algo.bayes currently can't handle Date entries. Computing reference values based on freq")
  }
  wrap.algo(sts,algo="algo.bayes",control=control)
}

#RKI wrapper
rki <- function(sts, control = list(range = range, b = 2, w = 4, actY = FALSE),...) {
  if (sts@epochAsDate) {
    warning("algo.rki currently can't handle Date entries. Computing reference values based on freq")
  }
  wrap.algo(sts,algo="algo.rki",control=control,...)
}

#outbreakP wrapper
outbreakP <- function(sts, control=list(range = range, k=100,
         ret=c("cases","value"),maxUpperboundCases=1e5),...) {
  wrap.algo(sts,algo="algo.outbreakP",control=control,...)
}


#HMM wrapper
hmm <- function(sts, control=list(range=NULL, noStates=2, trend=TRUE, noHarmonics=1,covEffectEqual=FALSE),...) {
  if (sts@epochAsDate) {
    warning("algo.hmm currently can't handle Date entries. Computing reference values based on freq")
  }
  wrap.algo(sts,algo="algo.hmm",control=control,...)
}


#Cusum wrapper
cusum <- function(sts,  control = list(range=range, k=1.04, h=2.26, m=NULL, trans="standard",alpha=NULL),...) {
  wrap.algo(sts,algo="algo.cusum",control=control,...)
}

#GLRpois wrapper
glrpois <- function(sts, control = list(range=range,c.ARL=5, S=1,
                           beta=NULL, Mtilde=1, M=-1, change="intercept",theta=NULL),...) {
  wrap.algo(sts,algo="algo.glrpois",control=control,...)
}

#GLRnb wrapper
glrnb <- function(sts, control = list(range=range,c.ARL=5, mu0=NULL, alpha=0,
                         Mtilde=1, M=-1, change="intercept",theta=NULL,dir=c("inc","dec"),
                         ret=c("cases","value")), ...) {
  wrap.algo(sts,algo="algo.glrnb",control=control,...)
}



#### this code definitely needs some more documentation -- wrap.algo atm is
# 100% without docu
#Rogerson wrapper
# theta0t now has to be a matrix
#library(surveillance)
#data("ha")
#rogerson(disProg2sts(ha),control=list(range=200:290,ARL0=100,s=1,theta0t=matrix(1,nrow=91,ncol=12)))

rogerson <- function(sts, control = list(range=range, theta0t=NULL,
                            ARL0=NULL, s=NULL, hValues=NULL,
                            distribution=c("poisson","binomial"),
                            nt=NULL, FIR=FALSE,limit=NULL, digits=1),...) {
  if (sts@epochAsDate) {
    warning("algo.rogerson currently can't handle Date entries. Computing reference values based on freq")
  }
  #Hook function to find right theta0t vector
  control.hook = function(k,control) {
    #Extract values relevant for the k'th component
    control$theta0t <- control$theta0t[,k]
    if (is.null(control[["nt",exact=TRUE]])) {
      control$nt <- sts@populationFrac[control$range,k]
    } else {
      if (!all.equal(sts@populationFrac[control$range,k],control$nt[,k])) {
        warning("Warning: nt slot of control specified, but specified population differs.")
      } else {
        control$nt <- control$nt[,k]
      }
    }
    #If no hValues given then compute them
    if (is.null(control[["hValues",exact=TRUE]])) {
#This code does not appear to work once n is big.
#      control$hValues <- hValues(theta0 = unique(control$theta0t), ARL0=control$ARL0, s=control$s , distr = control$distribution, n=mean(control$nt))$hValues
            control$hValues <- hValues(theta0 = unique(control$theta0t), ARL0=control$ARL0, s=control$s , distr = control$distribution)$hValues
    }
    return(control)
  }
  #WrapIt
  wrap.algo(sts,algo="algo.rogerson",control=control,control.hook=control.hook,...)
}


