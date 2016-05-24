###########################################################################
## A bunch of functions for point process simulations and analysis of
## these simulaions following the estimation of the intensity with gss. 
## Created: September 22 2009 by Christophe Pouzat
## Last modified: September 24 2009 by Christophe Pouzat
###########################################################################

maxIntensity <- function(object,
                         dfWithTime,
                         ...) {
############################################################################
## Simulates a point process (neuronal spike train) base on estimation
## of its intensity with gss. The simulation is performed with the thinning
## method.
## -------------------------------------------------------------------------
## Arguments:
##  object: A ssanova or ssanova0 object.
##  dfWithTime: A data frame with one variable named "time". The latter variable
##  is used to obtain the bin width with which the original spike train was
##  discretized.
##  ...: Additional arguments passed to optim which is called internally with
##  the "BFGS" method.
## --------------------------------------------------------------------------
## A "proposed" maximal intensity (in Hz) is returned.
###############################################################################  
  ## check that object inherists form ssanova or ssanova0
  if (!inherits(object, c("ssanova", "ssanova0"))) 
    stop("object should be a ssanova or a ssanova0 object.")
  ## check that the family used was binomial or poisson
  family4fit <- object[["call"]][["family"]]
  if (!(family4fit %in% c("binomial","poisson")))
    stop("The fit should have been done with the binomial or the poisson family.")
  mf <- object[["mf"]]
  allVN <- names(mf)[ names(mf) %in% object[["terms"]][["labels"]] ]
  nVar <- length(allVN)
  
  ## Get the bin width used to create the data frame to which the model
  ## was fitted.
  binWidth <- with(dfWithTime,diff(time)[1])

  nbPrev <- dim(mf)[1]
  df0 <- mf[nbPrev,allVN]
  ## Define IFct returning the intensity
  IFct <- switch(family4fit,
                 binomial=function(df=df0) {
                   pred <- exp(predict(object,df))
                   pred/(1+pred)/binWidth
                 },
                 poisson=function(df=df0) {
                   exp(predict(object,df))/binWidth
                 }
                 )

  targetFct <- function(p) {
    p <- exp(p);p <- p/(1+p)
    p <- matrix(p,nrow=1,dimnames=list(NULL,allVN))
    -predict(object,as.data.frame(p))
  }
  
  posMax <- optim(rep(0,length(allVN)),
                  targetFct,
                  method="BFGS",
                  ...)$par
  posMax <- exp(posMax);posMax <- posMax/(1+posMax)
  posMax <- matrix(posMax,nrow=1,dimnames=list(NULL,allVN))
  IFct(as.data.frame(posMax))*1.1
  
}

thinProcess <- function(object,
                        m2uFctList,
                        trueData,
                        formerSpikes,
                        intensityMax,
                        ...) {
############################################################################
## Simulates a point process (neuronal spike train) base on estimation
## of its intensity with gss. The simulation is performed with the thinning
## method.
## -------------------------------------------------------------------------
## Arguments:
##  object: A ssanova or ssanova0 object.
##  m2uFctList: A list of function. There should be as many functions as
##  there are "internal" variables in object. An internal variable is a variable
##  whose value is changed by the occurrence of a spike, like the elapsed
##  time since the last spike of the simulated neuron, the duration of a
##  former inter spike interval of a given lag, etc. The names of the
##  components (functions) of the list should be the names of the variables.
##  Each function should correspond to the map to uniform function used before
##  fitting the model.
##  trueData: A data frame containing the "true data" of the simulated epoch.
##  This is to ensure that "external" variables such as the elapsed time since
##  the last spike of a functionally coupled neuron are available.
##  formerSpikes: A vector of previous spike times. This is to make the
##  computation of former inter spike intervals possible in every case.
##  intensityMax: The value of the maximal intensity. If missing function
##  maxIntensity is called to estimate it.
##  ...: Additional argument for optim passed to maxIntensity if necessary.
## --------------------------------------------------------------------------
## A spike train object is returned.
###############################################################################

  ## check that object inherists form ssanova or ssanova0
  if (!inherits(object, c("ssanova", "ssanova0"))) 
    stop("object should be a ssanova or a ssanova0 object.")
  ## check that m2uFctList is a named list of functions
  if (is.null(names(m2uFctList)))
    stop("m2uFctList should be a NAMED list of functions.")
  if (!all(sapply(m2uFctList, function(f) inherits(f,"function"))))
    stop("m2uFctList should be a named list of FUNCTIONS.")
  family4fit <- object[["call"]][["family"]]
  if (!(family4fit %in% c("binomial","poisson")))
    stop("The fit should have been done with the binomial or the poisson family.")
  ##browser()
  mf <- object[["mf"]]
  allVN <- names(mf)[ names(mf) %in% object[["terms"]][["labels"]] ]
  nVar <- length(allVN)
  varVN <- allVN[allVN %in% names(m2uFctList)]

  ## Get the bin width used to create the data frame to which the model
  ## was fitted.
  binWidth <- with(trueData,diff(time)[1])

  nbPrev <- dim(mf)[1]
  df0 <- mf[nbPrev,allVN]
  ## Define IFct returning the intensity
  IFct <- switch(family4fit,
                 binomial=function(df=df0) {
                   pred <- exp(predict(object,df))
                   pred/(1+pred)/binWidth
                 },
                 poisson=function(df=df0) {
                   exp(predict(object,df))/binWidth
                 }
                 )
  
  if (missing(intensityMax)) {
    ## A maximal value for the intensity was not provided so
    ## one will be obtained by optimization
    intensityMax <- maxIntensity(object,trueData,...)
  } else {
    if (intensityMax<=0) stop("intensityMax must be > 0.")
  } ## End of conditional on missing(intensityMax)
  
  ## generate the realization of a uniform Poisson process with rate
  ## intensityMax
  from <- with(trueData,range(time))
  to <- from[2]
  from <- from[1]
  n2sim <- (to-from)*intensityMax
  n2sim <- ceiling(n2sim + 3*sqrt(n2sim))
  pProc <- from + cumsum(rexp(n2sim,intensityMax))
  while (max(pProc) < to) {
    pProc <- c(pProc, max(pProc) + cumsum(rexp(n2sim,intensityMax)))
  }
  pProc <- pProc[pProc <= to]
  ## Poisson process realization with rate intensityMax done
  
  st <- formerSpikes
  stLength <- length(st)
  for (poissonTime in pProc) {
    vVector <- sapply(m2uFctList,function(f) f(poissonTime,st))
    dfIdx <- (poissonTime-from) %/% binWidth + 1
    theDF <- trueData[dfIdx,]
    theDF[,varVN] <- vVector
    intensity <- IFct(theDF)
    if (intensity > intensityMax) {
      result <- list(IFct=IFct,
                     targetFct=targetFct,
                     posMax=posMax,
                     intensityMax=intensityMax)
      warning("intensityMax not large enough.")
      return(result)
    }
    if (runif(1,0,intensityMax) <= intensity) {
      ## A spike is generated
      st <- c(st,poissonTime)
    } ## End of conditional on runif(1,0,intensityMax) <= intensity
  } ## end of loop on poissonTime
  as.spikeTrain(st[-(1:stLength)])
}

mkSelf <- function(m2uSelf ## a map to uniform function
                           ## for the time to last spike
                   ) {
#############################################################################
## Create a function returning the mapped elapsed time since the last spike of
## the simulated neuron using the mapping to uniform which was used
## for the fitted part.
## --------------------------------------------------------------------------
## Argument:
##  m2uSelf: The map to uniform function used to transform the actual elapsed
##  time since the last spike values before fitting the model.
## --------------------------------------------------------------------------
## A function self(proposedtime,st) is returned.
#############################################################################
  force(m2uSelf)
  function(proposedtime,st) {
    m2uSelf(proposedtime-max(st))
  }
  
}


mkMappedI <- function(m2uI, ## a map to uniform function
                            ## for an interspike interval
                      lag = 1
                      ) {
#############################################################################
## Create a function returning a mapped former inter spike interval duration of
## the simulated neuron using the mapping to uniform which was used
## for the fitted part.
## --------------------------------------------------------------------------
## Argument:
##  m2uI: The map to uniform function used to transform the actual former
##  ISI durations before fitting the model.
##  lag: The considered lag (integer > 0).
## --------------------------------------------------------------------------
## A function(proposedtime,st) is returned.
#############################################################################  
  force(m2uI)
  force(lag)
  function(proposedtime,st) {
    n <- length(st)
    m2uI(diff(st[c(n-lag,n)]))
  }
  
}

mkSimFct <- function(object,
                     m2uFctList,
                     trueData,
                     formerSpikes,
                     intensityMax,
                     ...
                     ) {
#############################################################################
## Create a function returning a simulation by the thinning method of a
## spike train according to the model specified by object. The returned
## function is suitable for parallel use (eg, clusterApply) since it contains
## "everything needed" in its closure.
## --------------------------------------------------------------------------
## Argument:
##  object: A ssanova or ssanova0 object.
##  m2uFctList: A list of function. There should be as many functions as
##  there are "internal" variables in object. An internal variable is a variable
##  whose value is changed by the occurrence of a spike, like the elapsed
##  time since the last spike of the simulated neuron, the duration of a
##  former inter spike interval of a given lag, etc. The names of the
##  components (functions) of the list should be the names of the variables.
##  Each function should correspond to the map to uniform function used before
##  fitting the model.
##  trueData: A data frame containing the "true data" of the simulated epoch.
##  This is to ensure that "external" variables such as the elapsed time since
##  the last spike of a functionally coupled neuron are available.
##  formerSpikes: A vector of previous spike times. This is to make the
##  computation of former inter spike intervals possible in every case.
##  intensityMax: The value of the maximal intensity. If missing function
##  maxIntensity is called to estimate it.
##  ...: Additional argument for optim passed to maxIntensity if necessary.
## --------------------------------------------------------------------------
## A function simulating a spike train object is returned. The simulation is
## is done with function thinProcess. The returned function takes no argument.
## The maximal intensity required by the thinning method is stored in the
## closure of the returned function.
###############################################################################
  ## check that object inherists form ssanova or ssanova0
  if (!inherits(object, c("ssanova", "ssanova0"))) 
    stop("object should be a ssanova or a ssanova0 object.")
  ## check that the family used was binomial or poisson
  family4fit <- object[["call"]][["family"]]
  if (!(family4fit %in% c("binomial","poisson")))
    stop("The fit should have been done with the binomial or the poisson family.")
  ## check that m2uFctList is a named list of functions
  if (is.null(names(m2uFctList)))
    stop("m2uFctList should be a NAMED list of functions.")
  if (!all(sapply(m2uFctList, function(f) inherits(f,"function"))))
    stop("m2uFctList should be a named list of FUNCTIONS.")

  ## get names of functions in m2uFctList 
  m2uNames <- names(m2uFctList)
  ## create a list of environments, one per function in m2uFctList 
  envList <- lapply(m2uNames, function(fn) new.env())
  names(envList) <- m2uNames
  ## copy in each function's new environment what is in its present
  ## environment before setting each function's environment to
  ## the corresponding new environment
  for (fn in m2uNames) {
    ## get the names of the variables contained in
    ## each function's environment
    envVN <- ls(envir=environment(m2uFctList[[fn]]))
    lapply(envVN,
           function(vn) {
             exp <- parse(text=vn)
             value <- eval(exp,envir=environment(m2uFctList[[fn]]))
             assign(vn,value,envir=envList[[fn]])
           }
           )
    environment(m2uFctList[[fn]]) <- envList[[fn]]
  } ## End of loop on fn
  rm(envVN,fn)
  force(formerSpikes)
  force(trueData)
  
  if (missing(intensityMax)) {
    intensityMax <- maxIntensity(object,trueData,...)
  }
  function() thinProcess(object,
                         m2uFctList,
                         trueData,
                         formerSpikes,
                         intensityMax
                         )

}


mkPostSimAnalysis <- function(stList,
                              trainNumber=1,
                              timeWindow,
                              objects,
                              dfFct
                              ) {
#############################################################################
## Performs an automatic analysis, log predictive probability calculation and
## time transformation of simulated data. Several competing models can be used
## for the analysis.
## --------------------------------------------------------------------------
## Argument:
##  stList: The list of spikeTrain objects with one of the trains partly
##  simulated. A single (partly simulated) spikeTrain object can also be used.
##  trainNumber: An integer, the index of the modeled and simulated spike train
##  in stList.
##  timeWindow: A numeric vector of length 1 or 2. This argument specifies the
##  time domain over which the fits contained in argument objects was performed.
##  It is implicitly assumed that the (partial) simulation was performed outside
##  this time domain. When a vector of length 1 is used the fitting time domain
##  is taken as c(0,timeWindow).
##  objects: A list of ssanova or ssanova0 objects. Each element of the list
##  is a "model" with which analysis will be performed. A single ssanova or
##  ssanova0 object can also be used.
##  dfFct: a function whose argument is a the same as the first argument of
##  function mkGLMdf and which returns a data frame suitable for use as argument
##  newdata in predict.ssanova.
## --------------------------------------------------------------------------
## A list of lists is returned. Each list correspond to one of the models in
## argument objects. Each sub list has two components: lpp (the log predictive
## probability) and ttt (the time transformed train, a countingProcessSamplePath
## object).
###############################################################################  
  if (!inherits(stList,"list")) stList <- list(stList)
  earlyTrain <- unclass(stList[[trainNumber]])
  if (length(timeWindow) < 2) timeWindow <- c(0,timeWindow)
  earlyTrain <- earlyTrain[timeWindow[1] <= earlyTrain &
                           earlyTrain <= timeWindow[2]]
  if (!inherits(objects,"list")) objects <- list(objects)

  require(codetools)
  fNames <- findGlobals(dfFct,merge=FALSE)$functions
  fNamesG <- fNames[fNames %in% ls(globalenv())]
  if (length(fNamesG) > 0) {
    for (n in fNamesG) {
      assign(n,eval(parse(text=n)))
    }
  } ## end of conditional on length(fNamesG) > 0
  environment(dfFct) <- environment()
  function(st) { 
    fackedTrain <- as.spikeTrain(c(earlyTrain,unclass(st)))
    fackedSTList <- stList
    fackedSTList[[trainNumber]] <- fackedTrain
    DF <- dfFct(fackedSTList)
    lapply(objects,
           function(obj) {
             list(lpp = predictLogProb(obj,subset(DF,!(timeWindow[1] <= time &
                    time <= timeWindow[2]))),
                  ttt = obj %tt% subset(DF,!(timeWindow[1] <= time &
                    time <= timeWindow[2]))
                  )
           }
           )
  }
  
}
