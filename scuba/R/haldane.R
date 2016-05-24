#
# 	haldane.R
#
#	$Revision: 1.33 $	$Date: 2013/12/10 11:39:04 $
#
####################################################################
#
#  Calculations for Haldane type models
#

saturated.state <- function(model="D", depth=0, g=air) {
  if(is.character(model))
    model <- pickmodel(model)
  else
    stopifnot(inherits(model, "hm"))
  if(!capable(model, g))
    stop(paste("Model cannot deal with", as.character(g)))
  y <- summary(model)
  Pamb <- 1 + depth/10
  out <- list(N2= rep(Pamb * g$fN2, y$nc), He=rep(Pamb * g$fHe, y$nc))
  out <- out[y$species]
  out <- as.data.frame(out)
  row.names(out) <- summary(model)[["cnames"]]
  return(out)
}

conform <- function(state, model="D") {
  if(is.character(model))
    model <- pickmodel(model)
  else
  stopifnot(inherits(model, "hm"))
  m <- summary(model)
  mgas <- m$species
  s <-
    if(is.numeric(state) && is.vector(state)) data.frame(N2=state) else
    if(is.list(state)) as.data.frame(state) else
    if(is.data.frame(state)) state else
    if(is.matrix(state)) as.data.frame(state) else 
    stop("Unrecognised format for state")
  sgas <- colnames(s)
  if(is.null(sgas)) {
    if(ncol(s) == length(mgas))
      colnames(s) <- mgas
    else
      stop("gases in state do not match gases in model")
  } else {
    if(all(sgas %in% mgas) && all(mgas %in% sgas))
      s <- s[mgas]
    else
      stop("gases in state do not match gases in model")
  }
  if(nrow(s) != m$nc)
    stop("numbers of compartments in state do not match compartments in model")
  # attach names of compartments
  row.names(s) <- m[["cnames"]]
  return(s)
}

"haldane" <- 
function(d,
	 model=pickmodel("DSAT"),
	 prevstate=NULL,
         progressive=FALSE,
         relative=FALSE,
         deco=FALSE,
         derived=FALSE) 
{
  stopifnot(is.dive(d))

  if(is.character(model))
    model <- pickmodel(model)
  else
    stopifnot(inherits(model, "hm"))
  y <- summary(model)

  if(is.null(prevstate))
    prevstate <- saturated.state(model, 0)
  else
    prevstate <- conform(prevstate, model)

  if(relative) {
    if(!capable(model, d, "M0"))
      stop(paste("Cannot compute relative saturations:",
                 "the model does not include surfacing M-values",
                 "for some of the gases in this dive"))
    if(deco && !capable(model, d, "dM"))
      stop(paste("Cannot compute relative saturations for deco dive:",
                 "the model does not include M-value gradients (dM)",
                 "for some of the gases in this dive"))
  }

  species <- y$species
  
  # vectors: time
  times <- times.dive(d)
  depths <- depths.dive(d)
  durations <- diff(times)
  n <- length(times)

  # data frames: compartment x species
  state <- as.data.frame(prevstate)
  
  HalfT <- lapply(species,
                  function(x, ...) { param(species=x, what="HalfT", ...) },
                  model=model)
  names(HalfT) <- species
  HalfT <- data.frame(HalfT)
  ratecoefs <- log(2)/HalfT

  ncomp <- nrow(state)

  # data frame: time x species
  fGasList <- with(d$data, list(N2=fN2, He=fHe))
  fGas <- as.data.frame(fGasList[species])
  ngas <- ncol(fGas)

  # allocate space for result
  
  # ......... core calculation ...................................
  useC <- !identical(options("HaldaneC")[[1]], FALSE)
  if(useC) {
    # use C code
    pressures <- 1 + depths/10
    profile <- if(progressive) numeric(n * ncomp * ngas) else numeric(1)
    final <- numeric(n * ncomp * ngas)
    z <- .C("HaldaneCalc",
            ntimes   = as.integer(n),
            ntissues = as.integer(ncomp),
            ngases   = as.integer(ngas), 
            initial  = as.double(as.matrix(state)),
            pressures = as.double(pressures),
            durations = as.double(durations),
            gasfractions = as.double(as.matrix(fGas)),
            ratecoefs = as.double(as.matrix(ratecoefs)),
            progressive = as.integer(progressive),
            profile = as.double(profile),
            final = as.double(final)
            )
    state <- as.data.frame(matrix(z$final,
                                  nrow=nrow(state), ncol=ncol(state),
                                  dimnames=dimnames(state)))
    if(progressive) 
      profile <- array(z$profile, dim=c(n, dim(state)),
                       dimnames = append(list(NULL), dimnames(state)))
  } else {
    # interpreted code
    if(progressive) {
      profile <- array(0, dim=c(n, dim(state)))
      dimnames(profile) <- append(list(NULL), dimnames(state))
      profile[1,,] <- as.matrix(state)
    }
    for(i in 1:(n-1)) {
      d0 <- depths[i]
      d1 <- depths[i+1]
      tim <- durations[i]
      if(tim > 0) {
        fgi <- as.numeric(fGas[i,])
        k1 <- fgi * (d0/10 + 1)
        k2 <- fgi * ((d1-d0)/10) / tim
        k1 <- outer(rep(1, ncomp), k1, "*")
        k2 <- outer(rep(1, ncomp), k2, "*")
        v <- exp(- ratecoefs * tim)
        state <- state * v + ((k1 - k2/ratecoefs) * (1 - v) + k2 * tim)
      }
      if(progressive)
        profile[i+1,,] <- as.matrix(state)
    }
  }
  
  rownames(state) <- rownames(model)
  colnames(state) <- species

  if(derived) {
    # calculate derived quantities
    fGasAll <- as.matrix(as.data.frame(fGasList))
    if(!progressive) {
      # inspired gas partial pressure at end of dive
      ffGas <- matrix(fGasAll[n,], nrow=nrow(state), byrow=TRUE)
      dimnames(ffGas) <- dimnames(state)
      ppGas <- ffGas * (1+depths[n]/10)
      # washout at surface
      washout <- state - ppGas
      # decompression ceiling
      Pceiling <- deco.ceiling(profile, model, "pressure")
      Dceiling <- deco.ceiling(profile, model, "depth")
    } else {
      # inspired gas partial pressures (time x compartment x species)
      ffGas <- array(fGasAll, dim=c(n, rev(dim(state))))
      ffGas <- aperm(ffGas, c(1, 3, 2))
      dimnames(ffGas) <- dimnames(profile)
      ppGas <- ffGas * (1 + depths/10)
      # washout 
      washout <- profile - ppGas
      # decompression ceiling
      Pceiling <- deco.ceiling(profile, model, "pressure")
      Dceiling <- deco.ceiling(profile, model, "depth")
    }
    derivedstuff <- list(washout=washout,
                         Pceiling=Pceiling,
                         Dceiling=Dceiling)
  }

  if(!relative) {
    # absolute saturations
    result <- if(progressive) profile else state
  } else {
    # relative saturations
    fN2 <- d$data$fN2
    fHe <- d$data$fHe
    if(!progressive) {
      Mvals <- if(!deco) M0mix(model, fN2[n], fHe[n]) else
      Mmix(model, depths[n], fN2[n], fHe[n])    
      igstate <- apply(state, 1, sum)
      result <- relstate <- igstate/as.vector(Mvals)
    } else {
      Mvals <- if(!deco) M0mix(model,fN2,fHe) else Mmix(model,depths,fN2,fHe) 
      igprofile <- apply(profile, c(1,2), sum) 
      result <- relprofile <- igprofile/Mvals
    }
  }

  if(derived)
    attr(result, "derived") <- derivedstuff
  return(result)
}

###################################################################
#
#  Interactive display
#

showstates <- function(d, model="DSAT", relative=TRUE, deco=FALSE, ...) {
  dname <- deparse(substitute(d))
  stopifnot(is.dive(d))
  main <- list(...)$main
  if(is.null(main)) main <- dname

  # pick model by its name
  if(is.character(model))
    model <- pickmodel(model)
  else if(!inherits(model, "hm"))
    stop("model should be an object of class hm")

  if(!capable(model, d, "M0"))
    stop("The model does not include data for some of the gases in this dive")

  y <- summary(model)

  oldpar <- par(mfrow=c(1,2), ask=FALSE)
  frame()
  plot(d, ..., main=main)
  timepoints <- times.dive(d)
  maxtime <- max(timepoints)
  ntissues <- y$nc

  if(relative) {
    Mvals <- if(!deco) M0mix(model, d$data$fN2, d$data$fHe) else
                       Mmix(model, depths.dive(d), d$data$fN2, d$data$fHe)
  }
  

  # precompute time profile of tissue saturation
  cat("precomputing tissue saturations...")
  flush.console()
  halprofile <- haldane(d, model, progressive=TRUE)

  # sum over gas species N2 + He = 'inert gas'
  igprofile <- apply(halprofile, c(1,2), sum) 
                        
  if(relative) {
    # convert to relative saturation profile
    relprofile <- igprofile/Mvals
    # determine (approx) maximum relative saturation over entire dive
    ymax <- max(max(relprofile), 1)
  } else ymax <- max(max(igprofile))

  # precompute oxygen toxicity profile
  oxtoxprofile <- oxtox(d, progressive=TRUE)

  cat("done.\nplease click on the graph\n")
  flush.console()

  # event loop ..
  hal <- rep(NA, ntissues)

  while(length(xy <- locator(1)) > 0) {
    tim <- xy$x
    tim <- max(0, tim)
    tim <- min(tim, maxtime)

    ind  <- max(which(timepoints <= tim))
    ind1 <- min(which(timepoints >= tim))
    hal <- halprofile[ind,,] 
    ig  <- igprofile[ind,]
    if(relative) 
      rel <- as.numeric(relprofile[ind, ])
    otu <- oxtoxprofile[ind]

    depths <- depths.dive(d)
    if(ind < ind1) {
      # interpolate using correct equations
      t0 <- timepoints[ind]
      t1 <- timepoints[ind1]
      divebit <- chop.dive(dive.segment(d, ind), t0, tim)
      hal <- haldane(divebit, model, prevstate=hal)
      ig <- apply(hal, 1, sum)
      if(relative) {
        M. <- if(!deco) M0mix(model, d$data$fN2[ind], d$data$fHe[ind]) else
           Mmix(model, depths[ind], d$data$fN2[ind], d$data$fHe[ind]) 
        rel <- as.numeric(ig/M.)
      }
      otu <- otu + oxtox(divebit, warn=FALSE)
    }
    ylab <- if(!relative) "Tissue saturation (ata)" else
            if(!deco) "Relative saturation (for surfacing)" else
            "Relative saturation (at depth)" 
    pozzie <- barplot(if(relative) rel else ig,
                      xlab=paste("Tissues (", y$title, ")", sep=""),
                      ylab=ylab,
                      main=paste("Time=", round(tim), "min\n",
                        "Accumulated oxygen toxicity", round(otu, 1)),
                      ylim=range(c(0, 1.1 * ymax)),
                      names.arg=NULL)
    if(ntissues > 1) {
      mtext(side=1, at=pozzie[1], line=1, text="Fast")
      mtext(side=1, at=pozzie[ntissues], line=1, text="Slow")
    }
    abline(h=1, lty=3, col="red")
    plot(d, ..., main=dname)
    abline(v=tim, lty=2, col="green")
  }
  par(oldpar)
  if(is.matrix(hal))
    colnames(hal) <- y$species
  hal
}


###################################################################
#
#  NDL calculation
#

ndl <- function(depth, g=air, prevstate=NULL, model="DSAT") {
  stopifnot(is.gas(g))
  if(is.character(model))
    model <- pickmodel(model)
  else stopifnot(inherits(model, "hm"))
  if(!capable(model, g))
    stop(paste("Model cannot deal with", as.character(g)))
  if(is.null(prevstate))
    prevstate <- saturated.state(model)
  else
    prevstate <- conform(prevstate, model)
  result <- witch <- numeric(n <- length(depth))
  if(n == 0) return(result)
  species <- summary(model)$species
  pressure <- depth/10 + 1
  if(all(species == "N2")) {
    # N2 only
    fN2 <- g$fN2
    M0 <- param(model, "N2", "M0")
    ratecoefs <- log(2)/param(model, "N2", "HalfT")
    init <- prevstate$N2
    for(i in 1:n) {
      ppN2 <- fN2 * pressure[i]
      if(any(init > M0)) {
          # initial state of diver does not permit a dive
        result[i] <- 0
        witch[i] <- NA
      } else if(!any(bite <- (ppN2 > M0))) {
        # tissue tension will never exceed M-value in any compartment
        result[i] <- Inf
        witch[i] <- NA
      } else {
        stuff <- -log(((ppN2-M0)[bite])/((ppN2-init)[bite]))/ratecoefs[bite]
        result[i] <- min(stuff)
        witch[i] <- (seq(along=bite)[bite])[which.min(stuff)]
      }
    }
    attr(result, "controlling") <- witch
    return(result)
  } else {
    # N2 and He 
    fN2 <- g$fN2
    fHe <- g$fHe
    rateN2 <- log(2)/param(model, "N2", "HalfT")
    rateHe <- log(2)/param(model, "He", "HalfT")
    rateSLOWER <- pmin(rateN2, rateHe)
    ntissues <- length(rateN2)
    # compute surfacing M-value for gas
    M0 <- M0mix(model, fN2, fHe)
    # initial tissue saturations
    initN2 <- prevstate$N2
    initHe <- prevstate$He
    initIG <- initN2 + initHe
    #
    decay <- function(x, init, rate, asymp) {
      asymp + exp(-rate * x) * (init-asymp)
    }
    objective <- function(x, ppN2, initN2, rateN2, ppHe, initHe, rateHe, M0) {
      decay(x, initN2, rateN2, ppN2) + decay(x, initHe, rateHe, ppHe) - M0
    }
    epsilon <- .Machine$double.eps
    for(i in 1:n) {
      ppN2 <- fN2 * pressure[i]
      ppHe <- fHe * pressure[i]
      ppIG <- ppN2 + ppHe
      val <- numeric(ntissues)
      for(j in 1:ntissues) {
        if(initIG[j] > M0[j])
          # max tissue tension already exceeded
          val[j] <- 0
        else if(ppIG <= M0[j])
          # tissue tension will never exceed M-value in this compartment
          val[j] <-   Inf
        else {
          endpoint <- -log((ppIG-M0[j])/(ppIG-initIG[j]))/rateSLOWER[j]
          endvalue <- objective(endpoint, initN2=initN2[j],
                            rateN2=rateN2[j],
                            ppN2=ppN2,
                            initHe=initHe[j],
                            rateHe=rateHe[j],
                            ppHe=ppHe,
                            M0=M0[j])
          if(abs(endvalue) <= 2 * epsilon)
            val[j] <- endpoint
          else 
            val[j] <- uniroot(objective,
                              c(0, endpoint),
                              initN2=initN2[j],
                              rateN2=rateN2[j],
                              ppN2=ppN2,
                              initHe=initHe[j],
                              rateHe=rateHe[j],
                              ppHe=ppHe,
                              M0=M0[j])$root
        }
      }
      result[i] <- min(val)
      witch[i] <- which.min(val)
    }
    attr(result, "controlling") <- witch
    return(result)
  }
}


