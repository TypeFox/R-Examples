## Bayesian and bootstrap fitting of dose-response models

checkPrior <- function(prior){
  z <- 1
  for(z in 1:length(prior)){
    prvec <- prior[[z]]
    nam <- names(prior)[z]
    if(!all(is.numeric(prvec)))
      stop("non-numeric entry in prior")
    if(nam %in% c("norm", "t", "lnorm")){
      if(nam == "t"){
        if(length(prvec) != 3)
          stop("need vector of length 3 for ", nam, " prior")
        if(prvec[2] <= 0|prvec[3] <= 0)
          stop("2nd and 3rd entry needs to be positive for ", nam, " prior")
      } else {
        if(length(prvec) != 2)
          stop("need vector of length 2 for ", nam, " prior")
        if(prvec[2] <= 0)
          stop("2nd entry needs to be positive for ", nam, " prior")
      }
    } else {
      if(length(prvec) != 4)
        stop("need vector of length 4 for beta prior")
      if(min(prvec[3:4]) <= 0)
        stop("entry 3 and 4 need to be positive for beta prior")
      if(prvec[1] >= prvec[2])
        stop("entry 1 needs to be smaller than entry 2 for beta prior")
    }
  }
}

getPrBnds <- function(prior){
  prbnds <- matrix(ncol = 2, nrow = length(prior))
  for(z in 1:length(prior)){
    prvec <- prior[[z]]
    nam <- names(prior)[z]
    if(nam %in% c("norm", "t"))
      prbnds[z,] <- c(-Inf, Inf)
    if(nam == "lnorm")
      prbnds[z,] <- c(0, Inf)
    if(nam == "beta")
      prbnds[z,] <- c(prvec[1], prvec[2])
  }
  prbnds
}

projPrBnds <- function(par, lbnd, ubnd){
  ## project start parameter into bounds
  if(par > lbnd & par < ubnd){
    return(par)
  } else {
    rng <- ubnd-lbnd
    if(!is.finite(rng))
      rng <- 5
    if(par <= lbnd)
      return(lbnd+0.05*rng)
    if(par >= ubnd)
      return(ubnd-0.05*rng)
  }
}

bFitMod <- function(dose, resp, model, S, placAdj = FALSE,
                    type = c("Bayes", "bootstrap"),
                    start = NULL, prior = NULL, nSim = 1000,
                    MCMCcontrol = list(), control = NULL, bnds, 
                    addArgs = NULL){
  if(placAdj & model %in% c("linlog", "logistic"))
    stop("logistic and linlog models can only be fitted with placAdj")
  nD <- length(dose)
  if (length(resp) != nD) 
    stop("dose and resp need to be of the same size")
  dose <- as.numeric(dose)
  if (any(dose < -.Machine$double.eps)) 
    stop("dose values need to be non-negative")
  if (!is.numeric(dose)) 
    stop("dose variable needs to be numeric")
  resp <- as.numeric(resp)
  ## order dose and resp increasingly
  ord <- order(dose)
  dose <- dose[ord]
  resp <- resp[ord]
  if (nrow(S) != nD | ncol(S) != nD) 
    stop("S and dose have non-confirming size")
  if (missing(model)) 
    stop("need to specify the model that should be fitted")
  scal <- off <- nodes <- NULL
  if(model %in% c("linlog", "betaMod")){
    lst <- getAddArgs(addArgs, dose)
    if(model == "betaMod")
      scal <- lst$scal
    if(model == "linlog")
      off <- lst$off
  }
  if(model == "linInt")
    nodes <- dose
  
  ## model number
  builtIn <- c("linear", "linlog", "quadratic", "linInt", "emax",
               "logistic", "exponential", "sigEmax", "betaMod")
  modNr <- match(model, builtIn)
  if(is.na(modNr))
    stop("invalid model selected")
  ## number of parameters
  nPar <- as.integer(c(2, 2, 3, length(dose), 3, 4, 3, 4, 4)[modNr])

  type <- match.arg(type)
  if(type == "Bayes"){
    res <- bFitMod.Bayes(dose, resp, S, model, placAdj,
                         start, prior, nSim, MCMCcontrol,
                         off, scal, nPar, modNr)
    res <- matrix(res, nrow = nSim, ncol = nPar)
    if(placAdj & model != "linInt")
      res <- res[,-1, drop = FALSE]
  } else { ## bootstrap
    res <- bFitMod.bootstrap(dose, resp, S, model, placAdj,
                             nSim, control, bnds, off, scal,
                             nodes)
  }
  
  out <- list(samples = res)
  if(model != "linInt"){
    nams <- names(formals(model))[-1]
  } else {
    nams <- paste("d", dose, sep="")
  }
  if(modNr %in% c(2,9))
    nams <- nams[-length(nams)]
  if(placAdj & model != "linInt")
    nams <- nams[-1]
  colnames(out$samples) <- nams
  attr(out, "model") <- model
  lst <- list(dose, resp, S)
  doseNam <- as.list(match.call())$dose
  respNam <- as.list(match.call())$resp
  attr(out, "doseRespNam") <- as.character(c(doseNam, respNam))
  names(lst) <- c(doseNam, respNam, "S")
  attr(out, "data") <- lst
  attr(out, "type") <- type
  attr(out, "call") <- match.call()
  attr(out, "placAdj") <- placAdj
  attr(out, "prior") <- prior
  attr(out, "scal") <- scal
  attr(out, "off") <- off
  attr(out, "nodes") <- nodes
  class(out) <- "bFitMod"
  out
}

bFitMod.Bayes <- function(dose, resp, S, model, placAdj,
                          start, prior, nSim, MCMCcontrol,
                          off, scal, nPar, modNr){
  ## get defaults for MCMCcontrol
  ctrl <- list(thin = 1, w = NULL, adapt=TRUE)
  if (!is.null(MCMCcontrol)) {
    MCMCcontrol <- as.list(MCMCcontrol)
    ctrl[names(MCMCcontrol)] <- MCMCcontrol
  }
  ## check prior distribution
  if(is.null(prior))
    stop("need specification of prior in prior argument")
  prnr <- match(names(prior), c("norm", "t", "lnorm", "beta"))
  if(any(is.na(prnr)))
    stop("invalid prior selected")
  np <- nPar
  if(placAdj){
    if(model != "linInt"){ 
      np <- nPar - 1   
    } else {
      placAdj <- FALSE ## can proceed as if placAdj = FALSE
    }
  }
  if(length(prnr) != np)
    stop(length(prnr), " priors specified, need ", np," for selected model")
  checkPrior(prior)
  prBnds <- getPrBnds(prior)    
  
  ## add some checks here (scale > 0, a > b, alpha,beta>0)  
  prior <- as.double(do.call("c", prior))
  ## calculate starting value using fitMod if needed
  ## and width parameter for slice sampler
  if(is.null(start)|is.null(ctrl$w)){  
    mD <- max(dose)
    ll <- list(emax = c(0.1, 1.5) * mD, exponential = c(0.5, 1.5) * mD,
               logistic = matrix(c(0.1, 0.05, 1.5, 1/2) * mD, 2),
               sigEmax = matrix(c(0.1 * mD, 0.5, 1.5 * mD, 5), 2), 
               betaMod = matrix(c(0.2, 0.2, 4, 4), 2))
    gfit <- fitMod(dose, resp, S=S, model=model, type = "general",
                   bnds = ll[[model]],
                   placAdj = placAdj, addArgs=list(off = off, scal = scal))
    if(is.null(start)){
      start <- coef(gfit)
      for(i in 1:length(start)){
        start[i] <- projPrBnds(start[i], prBnds[i,1], prBnds[i,2])
      }
    } else {
      for(i in 1:length(start)){
        if((start[i] < prBnds[i,1]) | (start[i] > prBnds[i,2]))
          stop("specified start value not consistent with bounds on prior distribution")
      }
    }
    if(is.null(ctrl$w))
      ctrl$w <- rep(1.0, nPar)#sqrt(diag(vcov(gfit)))
  }
  if(np != length(start))
    stop("start of wrong length")
  if(placAdj){ # append 0
    if(model != "linInt")
      start <- c(0.0, start)
  }
  if(length(ctrl$w) != length(start))
    stop("w and start need to be of same size")
  ## add information for beta and linlog model
  if(model == "betaMod"){
    if(is.null(scal))
      stop("need scal parameter for betaMod")
    start <- c(start, as.double(scal))
  }
  if(model == "linlog"){
    if(is.null(off))
      stop("need off parameter for betaMod")
    start <- c(start, as.double(off))
  }
  ## preliminary formatting to send information to C
  start <- as.double(start)
  inS <- solve(S)
  if(inherits(inS, "try-error")) 
    stop("specified S is not invertible")
  clinS <- as.double(chol(inS))
  ## ensure that parameters are of right class
  nSimTot <- as.integer(nSim*ctrl$thin);thin <- as.integer(ctrl$thin)
  out <- double(floor(nSimTot/thin)*nPar)  
  resp <- as.double(resp);dose <- as.double(dose)
  modNr <- as.integer(modNr);clinS <- as.double(clinS)
  nD <- as.integer(length(dose));w <- as.double(ctrl$w)
  noint <- as.integer(placAdj)
  ## call c code
  if(ctrl$adapt){
    res <- .C("sample", as.integer(500), as.integer(1), out=double(500*nPar),
              start, noint, w, dose, modNr, nPar, double(length(dose)),
              resp, clinS, nD, prior, prnr, double(nPar), double(nPar))
    res <- matrix(res$out, nrow = 500, ncol = nPar)
    w <- apply(res, 2, function(x) IQR(x)/1.3)
  }
  res <- .C("sample", nSimTot, thin, out=out, start, noint, w,
            dose, modNr, nPar, double(length(dose)), resp, clinS,
            nD, prior, prnr, double(nPar), double(nPar))
  res$out
}

bFitMod.bootstrap <- function(dose, resp, S, model, placAdj,
                              nSim, control, bnds, off, scal,
                              nodes){
  if(model %in% c("emax", "exponential", "betaMod", "logistic", "sigEmax")){
    if(missing(bnds)){
      message("Message: Need bounds in \"bnds\" for nonlinear models, using default bounds from \"defBnds\".")
      bnds <- defBnds(max(dose))[[model]]
    }
  }

  ## same arguments as in gFitDRModel function
  sims <- rmvnorm(nSim, resp, S)
  func <- function(x){
    fit <- fitMod.raw(dose, x, S=S, model=model, type="general",
                      placAdj=placAdj, bnds=bnds, control=control,
                      off=off, scal=scal, nodes=nodes,
                      covarsUsed = FALSE, df = Inf,
                      doseNam = "dose", respNam = "resp")
    coef(fit)
  }
  out <- apply(sims, 1, func)
  if(is.matrix(out)){
    return(t(out))
  } else {
    return(matrix(out, nrow = nSim, ncol = 1))
  }
}

## to do write print, predict and summary method
ess.mcmc <- function(series, lag.max = NULL){
  ## initial monotone sequence estimate of effective sample size
  ## Geyer, 1992, Statistical Science, idea:
  ## sum of even and un-even autocorrelations (gamma)
  ## needs to be positive and monotone decreasing
  N <- length(series)
  if (length(unique(series)) == 1) 
    return(NA)
  if (is.null(lag.max)) 
    lag.max <- 10 * log10(N)
  ac <- acf(series, plot = FALSE, lag.max = lag.max)$acf[2:(lag.max + 
                                    1), , 1]
  gam <- ac[-length(ac)]+ac[-1]
  dgam <- -diff(gam)
  if (gam[1] < 0) 
    return(N)
  m1 <- m2 <- lag.max
  ind1 <- gam < 0
  if (any(ind1)) 
    m1 <- min(which(ind1))
  ind2 <- dgam < 0
  if (any(ind2)) 
    m2 <- min(which(ind2))
  ind <- min(2 * min(m1, m2) + 1, lag.max)
  N/(1 + 2 * sum(ac[1:ind]))
}

print.bFitMod <- function(x, digits = 3, ...){
  ## print brief summary of MCMC samples
  doseNam <- attr(x, "doseRespNam")[1]
  respNam <- attr(x, "doseRespNam")[2]
  resp <- attr(x, "data")[[respNam]]
  names(resp) <- attr(x, "data")[[doseNam]]
  cat("Dose Response Model\n\n")
  cat(paste("Model:", attr(x, "model")), "\n\n")
  if(attr(x, "type") == "Bayes"){
    cat("Summary of posterior draws\n")
    func <- function(x){
      c(mean=mean(x), sdev=sd(x),
        quantile(x, c(0.025, 0.25, 0.5, 0.75, 0.975)),
        n.eff=ess.mcmc(x))
    }
    print(t(apply(x$samples, 2, func)), digits=digits)
  } else {
    cat("Summary of bootstrap draws\n")
    func <- function(x){
      c(mean=mean(x), sdev=sd(x),
        quantile(x, c(0.025, 0.25, 0.5, 0.75, 0.975)))
    }
    print(t(apply(x$samples, 2, func)), digits=digits)
  }
  cat("\nFitted to:\n")
  print(signif(resp, digits+2))
}

predict.bFitMod <- function(object, predType = c("full-model", "effect-curve"),
                            summaryFct = function(x) quantile(x, probs = c(0.025, 0.25, 0.5, 0.75, 0.975)),
                            doseSeq = NULL, lenSeq = 101, ...){
  predType <- match.arg(predType)
  doseNam <- attr(object, "doseRespNam")[1]
  if (is.null(doseSeq)) {
    doseSeq <- seq(0, max(attr(object, "data")[[doseNam]]), length = lenSeq)
  }
  model <- attr(object, "model")
  scal <- attr(object, "scal")
  off <- attr(object, "off")
  placAdj <- attr(object, "placAdj")
  if(placAdj){
    nodes <- c(0,attr(object, "data")[[doseNam]])
  } else {
    nodes <- attr(object, "data")[[doseNam]]
  }
  out <- predSamples(samples = object$samples, doseSeq = doseSeq,
                     placAdj = placAdj, model = model, scal = scal,
                     off = off, nodes = nodes)

  if(predType == "effect-curve"){
    out <- out - out[,1]
  }
  if(!is.null(summaryFct)){
    out0 <- apply(out, 2, summaryFct)
    out <- matrix(out0, ncol = ncol(out))
  }
  colnames(out) <- doseSeq
  out
}

predSamples <- function(samples, placAdjfullPars = FALSE, doseSeq, placAdj, model,
                        scal, off, nodes){
  ## placAdjfullPars argument only of interest if placAdj = TRUE
  ## it determines whether the e0 parameter is included as a row in the
  ## samples argument or not
  if(model != "betaMod")
    scal <- NULL
  if(model != "linlog")
    off <- NULL
  if(placAdj){
    if(placAdjfullPars){
      if(model != "linInt"){
        func <- function(x){
          pred <- do.call(model, c(list(doseSeq), as.list(c(x, scal, off))))
          pred0 <- do.call(model, c(list(0), as.list(c(x, scal, off))))
          pred-pred0
        }
      } else {
        func <- function(x){
          pred <- do.call(model, c(list(doseSeq), as.list(list(x, nodes))))
          pred0 <- do.call(model, c(list(0), as.list(list(x, nodes))))
          pred-pred0
        }
      }
    } else {
      if(model != "linInt"){
        func <- function(x)
          do.call(model, c(list(doseSeq), as.list(c(c(0,x), scal, off))))
      } else {
        func <- function(x)
          do.call(model, c(list(doseSeq), as.list(list(c(0,x), nodes))))
      }
    }
  } else {
    if(model != "linInt"){
      func <- function(x)
        do.call(model, c(list(doseSeq), as.list(c(x, scal, off))))
    } else {
      func <- function(x)
        do.call(model, c(list(doseSeq), as.list(list(x, nodes))))
    }
  }
  out <- t(apply(samples, 1, func))
}

plot.bFitMod <- function (x, plotType = c("dr-curve", "effect-curve"),
                          quant = c(0.025, 0.5, 0.975), 
                          plotData = c("means", "meansCI", "none"),
                          level = 0.95, lenDose = 201, ...){
  addArgs <- list(...)
  plotType <- match.arg(plotType)
  doseNam <- attr(x, "doseRespNam")[1]
  respNam <- attr(x, "doseRespNam")[2]
  dose <- attr(x, "data")[[doseNam]]
  resp <- attr(x, "data")[[respNam]]
  doseSeq <- seq(0, max(dose), length = lenDose)
  plotData <- match.arg(plotData)
  placAdj <- attr(x, "placAdj")
  sumFct <- function(x){
    quantile(x, probs = quant)
  }
  if (placAdj) 
    plotType <- "effect-curve"
  if (plotType == "effect-curve") {
    pred <- predict(x, predType = plotType, summaryFct = sumFct,
                    doseSeq = doseSeq)
    main <- "Effect Curve"
    if (placAdj) {
      if (plotData == "meansCI") {
        sdev <- sqrt(diag(attr(x, "data")$S))
        q <- qnorm(1 - (1 - level)/2)
        LBm <- UBm <- numeric(length(dose))
        for (i in 1:length(dose)) {
          LBm[i] <- resp[i] - q * sdev[i]
          UBm[i] <- resp[i] + q * sdev[i]
        }
      }
      else {
        LBm <- UBm <- NULL
      }
    }
    else {
      LBm <- UBm <- NULL
    }
  }
  if (plotType == "dr-curve") {
    pred <- predict(x, predType = "full-model", summaryFct = sumFct,
                    doseSeq = doseSeq)
    main <- "Dose-Response Curve\n"
    if (plotData == "meansCI") {
      sdev <- sqrt(diag(attr(x, "data")$S))
      q <- qnorm(1 - (1 - level)/2)
      LBm <- UBm <- numeric(length(dose))
      for (i in 1:length(dose)) {
        LBm[i] <- resp[i] - q * sdev[i]
        UBm[i] <- resp[i] + q * sdev[i]
      }
    }
    else {
      LBm <- UBm <- NULL
    }
  }
  rng <- range(c(pred, resp, LBm, UBm))
  dff <- diff(rng)
  ylim <- c(rng[1] - 0.02 * dff, rng[2] + 0.02 * dff)
  callList <- list(doseSeq, t(pred), type = "l", xlab = doseNam, 
                   ylim = ylim, ylab = respNam, main = main,
                   lty=1, col=1)
  callList[names(addArgs)] <- addArgs
  do.call("matplot", callList)
  
  if (plotType == "dr-curve" | placAdj) {
    if (plotData == "means") 
      points(dose, resp, pch = 19, cex = 0.75)
    if (plotData == "meansCI") {
      points(dose, resp, pch = 19, cex = 0.75)
      for (i in 1:length(dose)) {
        lines(c(dose[i], dose[i]), c(LBm[i], UBm[i]), 
              lty = 2)
      }
    }
  }
  res <- list(doseSeq = doseSeq)
  attr(res, "level") <- level
  attr(res, "ylim") <- ylim
  res$mean <- pred
  invisible(res)
}

coef.bFitMod <- function (object, ...){
  object$samples
}
