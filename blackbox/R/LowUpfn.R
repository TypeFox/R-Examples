LowUpfn <- function(lower, upper, boundstype="logical", scale=.blackbox.data$options$samplingScale,verbose=FALSE) {
  ## assumes lower, upper are in samplingSpace (necess to apply the constraints);
  ## to constrain extrapolated points
  ## expands (allowing obs/2 minimal values [in canonical scale, cannot be 0 if we relog afterwards] and Inf maximal values),
  ## with exceptions for more constrained parameters
  DemographicModel <- blackbox.getOption("DemographicModel")
  if ("IBD" %in% DemographicModel) D2bool <- ("2D" %in% DemographicModel)
  if (boundstype=="logical") {
    maxima <- list(g=1, pGSM=1, Q1=0.999999999,twoNm=Inf, twoNmu=Inf, twoNfoundermu=Inf, twoNancmu=Inf, `T`=Inf, D=Inf,
                   condS2=Inf)
    minima <- list(g=0, pGSM=0, Q1=0.000000001,twoNm=0, twoNmu=0, twoNfoundermu=0, twoNancmu=0, `T`=0, D=0,
                   condS2=condaxialS2fromg(0, D2bool=D2bool) )
  } else if (boundstype=="numerical") {
    maxima <- list(pGSM=0.99,Q1=0.999999999)
    minima <- list(g=0,pGSM=0.01,Q1=0.000000001)
    if ("IBD" %in% DemographicModel) {
      locnames <- unlist((lapply(blackbox.getOption("LRTlist"), function(l) {l$LRTnames})))
      if (length(intersect(c("condS2", "latt2Ns2"), c(blackbox.getOption("FONKgNames"), locnames)))>0) {
        maxima$g <- 0.999999 ## g=1 = infinite Nb ##FR->FR g=1 is not a pb for likelhood computation but if condS2 or Nb is used in kriging or LRT
      } else maxima$g <- 1
      maxima$condS2 <- condaxialS2fromg(maxima$g, D2bool=D2bool)
      minima$condS2 <- condaxialS2fromg(minima$g, D2bool=D2bool)
      positivemin <- 1e-06
    } else {positivemin <- 1e-09}
    minima <- c(minima, list(twoNm=positivemin,
                             twoNmu=positivemin,
                             twoNfoundermu=positivemin,
                             twoNancmu=positivemin,
                            `T`=positivemin,
                             D=positivemin))
    if(blackbox.getOption("LikStatistic") %innc% c("IS","ISstrict")) { ## FR->FR test a sortir de la fn?
      ## large 2Nm values *are* meaningful; but not feasible with IS
      maxima <- c(maxima, list(twoNm=1000, twoNmu=1000, twoNfoundermu=1000, twoNancmu=10000))
    } else {maxima <- c(maxima, list(twoNm=Inf, twoNmu=Inf, twoNfoundermu=Inf, twoNancmu=Inf))}
    # ajustements par option utilisateur
    allmaxnames <- names(maxima)
    allminnames <- names(minima)
    maxima <- lapply(allmaxnames, function(st) min(maxima[[st]], blackbox.getOption("ParMaxima")[[st]]))
    minima <- lapply(allminnames, function(st) max(minima[[st]], blackbox.getOption("ParMinima")[[st]]))
    names(minima) <- allminnames
    names(maxima) <- allmaxnames
  } else stop.redef("(!) From LowUpfn: unknown boundstype.")
  localLow <- lower;localUp <- upper
  ## unlogs
  logs <- tolower(scale)=="logscale" ## vector of T/F ## names of scale should be those of lower => default scale must be samplingScale (modified 20/06/13)
  logs <- logs[names(lower)] ## names determined by BoundsDefs
  localLow[logs] <- exp(localLow[logs])
  localUp[logs] <- exp(localUp[logs])
  ##shifts bounds
  for (st in names(lower)) {
    ## localLow ## FR->FR. Bounds for checking points after expansion by extrapol. So the can be rather lax. BUT
    ## FR->FR probably a good idea to make them depend on testPoint value, and to rename them
    if(st== "condS2") {
      if ("2D" %in% DemographicModel) {machin <- 1/2} else {machin <- 1}
      localLow[st] <- machin+(localLow[st]-machin)/2
    } else if (st %in% c("pGSM", "Q1","twoNm", "twoNmu", "twoNancmu", "twoNfoundermu", "T", "D")) {
      localLow[st] <- max(localLow[st]/50, minima[[st]]) ## cautious approach to zero.
    } else localLow[st] <- max(localLow[st]/2, minima[[st]])
    # # localUp ## basic idea is twofold expansion (in canonical scale)
    if (st == "g") {
      localUp[st] <- min((1+localUp[st])/2, maxima[[st]])
    } else if (st %in% c("pGSM", "Q1", "twoNm", "twoNmu", "twoNancmu", "twoNfoundermu", "condS2", "T", "D")) {
      localUp[st] <- min(2*localUp[st], maxima[[st]]) ## this may be large if high Nm value returned by a LRT , cf comments of purging locLRTlist
    } else {localUp[st] <- 2*localUp[st]}
  }
  userNLow <- localLow;names(userNLow) <- sapply(names(userNLow), formatName, format="ASCII")
  userNUp <- localUp;names(userNUp) <- sapply(names(userNUp), formatName, format="ASCII")
  if (verbose) {
    message.redef("Lower and upper allowed values of points to be generated: ")
    message.redef(rbind(prettynum(userNLow), prettynum(userNUp)))
  }
  ##relogs
  localLow[logs] <- log(localLow[logs])
  localUp[logs] <- log(localUp[logs])
  return(list(localLow=localLow, localUp=localUp))
} ## end def LowUpFn
