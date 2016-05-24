## various functions for assessing the operating characteristics of a design
## for model-based estimation of dose-response functions

## calculates the variance of the estimated curve
getPredVar <- function(model, cf, V,  pDose, off, scal){
  gr <- gradCalc(model, cf, pDose, off, scal)
  gr0 <- gradCalc(model, cf, 0, off, scal)
  grd <- sweep(gr, 2, gr0)
  out <- apply(grd, 1, function(x){
    as.numeric(t(x)%*%V%*%x)
  })
  out
}

## calculates the variance of the EDp estimate
getEDVar <- function(model, cf, V, scale = c("unrestricted", "logit"),
                     p, maxD, off, scal, nodes){
  grd <- calcEDgrad(model, cf, maxD, p, off, scal, nodes)
  if(scale == "logit"){
    tmp <- calcED(model, cf, p, maxD, "continuous",
                                off=off, scal=scal, nodes=nodes)
    grd <- grd*(-maxD/(tmp*(tmp-maxD)))
  }
  grd <- as.numeric(grd)
  return(as.numeric(t(grd)%*%V%*%grd))
}

## calculates the variance of the TD estimate
getTDVar <- function(model, cf, V, scale = c("unrestricted", "log"),
                     Delta, direction = c("increasing", "decreasing"), off, scal, nodes){
  tmp <- calcTD(model, cf, Delta, "continuous", 
                              direction, off=off, scal=scal,
                              nodes = nodes)
  grd <- calcTDgrad(model, cf, Delta, direction, off, scal, nodes)
  if(scale == "log")
    grd <- grd/tmp
  grd <- as.numeric(grd)
  return(as.numeric(t(grd)%*%V%*%grd))
}

## calculates approximate covariance matrix for parameter estimates
aprCov <- function(doses, model, cf, S, off, scal){
  F <- gradCalc(model, cf, doses, off, scal)
  V <- try(solve(t(F)%*%solve(S)%*%F))
  if(inherits(V, "try-error")){
    warning("Could not calculate covariance matrix; Fisher information singular.")
    return(NA)
  }
  V
}

planMod <- function(model, altModels, n, sigma, S, doses,
                    asyApprox = TRUE, simulation = FALSE,
                    alpha = 0.025, tau = 0,
                    p = 0.5, pLB = 0.25, pUB = 0.75,
                    nSim = 100, cores = 1,  showSimProgress = TRUE,
                    bnds, addArgs = NULL){
  if(any(is.element(model, "linInt")))
    stop("planMod works for all built-in models but not linInt")
  if(length(model) > 1 & asyApprox){
    stop("\"asyApprox\" needs to be FALSE for multiple models")
  } 
  ## off and scal
  off <- scal <- NULL
  if(any(is.element(model, c("linlog", "betaMod")))) {
    lst <- getAddArgs(addArgs, sort(unique(doses)))
    if ("betaMod" %in% model) 
      scal <- lst$scal
    if ("linlog" %in% model) 
      off <- lst$off
  }
  if(missing(doses))
    doses <- attr(altModels, "doses")
  ## calculate mean response at doses
  muMat <- getResp(altModels, doses)

  nD <- length(doses)
  if(missing(S)){
    if(missing(n) | missing(sigma))
      stop("either S or n and sigma need to be specified")
    if (length(n) == 1) 
      n <- rep(n, nD)
    if (length(n) != nD) 
      stop("\"n\" and \"doses\" need to be of same length")
    S <- sigma^2 * diag(1/n)
  }
  
  ## calculate parameters, gradients and results for the asymptotic approximation
  if(missing(bnds)) {
    if(any(!is.element(model, c("linear", "linlog", "quadratic")))){
      message("Message: Need bounds in \"bnds\" for nonlinear models, using default bounds from \"defBnds\".")
      bnds <- defBnds(max(doses))
    }
  }
  nams <- colnames(muMat)
  covMat <- list()
  approx <- matrix(nrow = ncol(muMat), ncol = 3)
  maxdose <- apply(abs(muMat-muMat[1,]), 2, function(x) doses[which.max(x)])
  EDs <- ED(altModels, p)
  EDsUB <- ED(altModels, pUB)
  EDsLB <- ED(altModels, pLB)

  if(!asyApprox & !simulation)
    stop("Need to select either \"asyApprox = TRUE\" or \"simulation = TRUE\"")
  
  if(asyApprox){
    npar <- switch(model,
                   linInt = length(doses),
                   nPars(model))
    bestPar <- matrix(nrow = ncol(muMat), ncol = npar) ## best fit by model to models in altModels
    for(i in 1:ncol(muMat)){
      ## if other model-class approximate best fit
      nam <- gsub("[0-9]", "", nams[i]) # model name (number removed)
      if(nam == model){
        pars <- attr(muMat, "parList")[[i]]
        if(is.element(model, c("betaMod", "linlog")))
          bestPar[i,] <- pars[-length(pars)]
        else
          bestPar[i,] <- pars
        bias <- 0
      } else { ## find the best fit 
        fit <- fitMod(doses, muMat[,i], model=model, S=S,
                      bnds = bnds[[model]], type="general")
        bias <- predict(fit, predType = "effect-curve" , doseSeq = doses[-1])-(muMat[-1,i]-muMat[1,i])
        bestPar[i,] <- coef(fit)
      }
      ## now calculate approximate covariance matrix
      covMat[[i]] <- aprCov(doses, model, bestPar[i,], S, off, scal)
      if(!is.matrix(covMat[[i]])){
        approx[i,] <- NA
      } else {
        ## root-mse
        paVar <- getPredVar(model, bestPar[i,], covMat[[i]], 
                            pDose=doses[-1], scal=scal, off=off)
        approx[i,1] <- sqrt(mean(paVar+bias^2))
        ## Pr(eff_maxdose > 0)
        ind <- which(doses[-1] == maxdose[i])
        paVar <- paVar[ind]
        call <- c(list(c(0,maxdose[i])), as.list(c(bestPar[i,], scal, off)))
        pa <- abs(diff(do.call(model, call)))
        LBmn <- qnorm(alpha, pa, sqrt(paVar))
        approx[i,2] <- pnorm(tau, LBmn, sqrt(paVar), lower.tail = FALSE)
        ## Pr(eff_ED50)
        edvar <- getEDVar(model, bestPar[i,], covMat[[i]], "unrestricted", p,
                          maxdose[i], off=off, scal=scal)
        ed <- calcED(model, bestPar[i,], p, maxdose[i], "continuous",
                                   off=off, scal=scal)
        edsd <- sqrt(edvar)
        approx[i,3] <- pnorm(EDsUB[i], ed, edsd) - pnorm(EDsLB[i], ed, edsd)
      }
    }
    colnames(approx) <- c("dRMSE", "Pow(maxDose)", "P(EDp)")
    rownames(approx) <-   rownames(bestPar) <- nams
    colnames(bestPar) <- rownames(covMat[[1]])
    attr(approx, "bestPar") <- bestPar
    attr(approx, "covMat") <- covMat
  }
  
  if(simulation){
    cat("Running simulations\n")
    requireNamespace("parallel", quietly = TRUE)
    sim <- parallel::mclapply(1:ncol(muMat), function(i){
      if(showSimProgress){
        if(cores == 1){
          cat(sprintf("Scenario %d/%d\n", i, ncol(muMat)))
          pb <- txtProgressBar(style=3, char="*")
        } else {
          cat(sprintf("Scenario %d/%d started\n", i, ncol(muMat)))
        }
      }
      dat <- rmvnorm(nSim, mean = muMat[,i], sigma = S)
      sims <- numeric(3)
      mse <- LBmn <- edpred <- resp <- numeric(nSim)
      coefs <- vector("list", length = nSim)
      modelSel <- character(nSim)
      for(j in 1:nSim){
        if(showSimProgress & cores == 1)
          setTxtProgressBar(pb, j/nSim)
        fit <- vector("list", length = length(model))
        k <- 1
        for(namMod in model){
          fit[[k]] <- fitMod(dose=doses, dat[j,], model=namMod,
                        S=S, type="general", bnds=bnds[[namMod]])
          k <- k+1
          ## ## this would be faster
          ## fit <- fitMod.raw(doses, dat[j,], model=model,
          ##                                 off=off, scal=scal, nodes=NULL,
          ##                                 S=S, type="general", bnds=bnds,
          ##                                 covarsUsed = FALSE, df = Inf,
          ##                                 control = NULL, 
          ##                                 doseNam = "dose", respNam = "resp")
        }
        aics <- sapply(fit, gAIC)
        fit <- fit[[which.min(aics)]]
        coefs[[j]] <- coef(fit)
        modelSel[j] <- attr(fit, "model")
        
        ## root-MSE of plac-adj dr at doses
        respDoses <- predict(fit, predType = "effect-curve", doseSeq = doses[-1])
        call <- c(list(doses), as.list(c(coef(fit), scal, off)))
        trm <- muMat[-1,i] - muMat[1,i]
        mse[j] <- mean((respDoses-trm)^2)
        ## Pr(LB_maxdose > tau) > 1-alpha
        respMaxD <- predict(fit, predType = "effect-curve", doseSeq = maxdose[i], se.fit=TRUE)
        if(is.na(respMaxD$se.fit)){
          LBmn[j] <- NA
        } else {
          LBmn[j] <- qnorm(alpha, abs(respMaxD$fit), respMaxD$se.fit)
        }
        resp[j] <- respMaxD$fit
        ## ED estimation
        edpred[j] <- ED(fit, p=p)
      }
      ind <- is.na(LBmn)
      NAind <- sum(ind)
      LBmn[ind] <- qnorm(alpha, abs(resp[ind]), sd(resp, na.rm=TRUE))
      sims[1] <- sqrt(mean(mse))
      sims[2] <- mean(LBmn > tau)
      sims[3] <- mean(edpred > EDsLB[i] & edpred < EDsUB[i])
      attr(sims, "NAind") <- NAind
      attr(sims, "coefs") <- coefs
      attr(sims, "model") <- modelSel 
      if(showSimProgress){
        if(cores == 1){
          close(pb)
        } else {
          cat(sprintf("Scenario %d/%d finished\n", i, ncol(muMat)))
        }
      }
      sims
    }, mc.cores=cores)
    NAind <- sapply(sim, function(x) attr(x, "NAind"))
    coefs <- lapply(sim, function(x) attr(x, "coefs"))
    modelSel <- sapply(sim, function(x) attr(x, "model"))
    names(NAind) <- colnames(modelSel) <- names(coefs) <- nams
    rownames(modelSel) <- 1:nSim
    sim <- do.call("rbind", sim)
    colnames(sim) <- c("dRMSE", "Pow(maxDose)", "P(EDp)")
    rownames(sim) <- nams
    attr(sim, "NAind") <- NAind
    attr(sim, "coefs") <- coefs
    attr(sim, "modelSel") <- modelSel
  }

  out <- list(approx = NULL, sim = NULL)
  if(asyApprox)
    out$approx <- approx
  if(simulation){
    out$sim <- sim
    attr(out$sim, "nSim") <- nSim
  }
  attr(out, "model") <- model
  attr(out, "altModels") <- altModels
  attr(out, "doses") <- doses
  attr(out, "off") <- off
  attr(out, "scal") <- scal
  attr(out, "S") <- S
  class(out) <- "planMod"
  out
}

tableMatch <- function(x, match){
  ## like "table", but also returns categories with 0 counts
  out <- numeric(length(match))
  for(i in 1:length(match)){
    out[i] <- sum(x == match[i], na.rm=TRUE)
  }
  names(out) <- match
  out
}


print.planMod <- function(x, digits = 3,...){
  model <- attr(x, "model")
  multiMod <- length(model) > 1
  str <- ifelse(multiMod, "s", "")
  cat(sprintf("Fitted Model%s: %s\n\n", str, paste(model, collapse=" ")))
  if(!is.null(x$approx)){
    attr(x$approx, "bestPar") <- NULL
    attr(x$approx, "NAind") <- NULL
    attr(x$approx, "covMat") <- NULL
    cat("Asymptotic Approximations\n")
    print(signif(x$approx, digits))
    cat("\n")
  }
  if(!is.null(x$sim)){
    pltsim <- x$sim
    attr(pltsim, "NAind") <- NULL
    attr(pltsim, "coefs") <- NULL
    attr(pltsim, "modelSel") <- NULL
    attr(pltsim, "nSim") <- NULL
    cat(sprintf("Simulation Results (nSim = %i)\n", attr(x$sim, "nSim")))
    print(signif(pltsim, digits))
    if(multiMod){
      cat("\nSelected models\n")
      res <- apply(attr(x$sim, "modelSel"), 2, tableMatch, match = model)
      print(signif(t(res)/colSums(res), digits))
    }
  }
}

summary.planMod <- function(object, digits = 3, len = 101,
                            Delta=NULL, direction = c("increasing", "decreasing"),
                            p=NULL, dLB = 0.05, dUB = 0.95, ...){
  class(object) <- "summary.planMod"
  print(object, digits, len, Delta, direction,
        p, dLB, dUB, ...)
}

print.summary.planMod <- function(x, digits = 3, len = 101,
                                  Delta=NULL, direction = c("increasing", "decreasing"),
                                  p=NULL, dLB = 0.05, dUB = 0.95, ...){
  ## provide more information than print method
  modelSel <- attr(x$sim, "modelSel") 
  model <- attr(x, "model")
  coefs <- attr(x$sim, "coefs")
  altModels <- attr(x, "altModels")
  doses <- attr(x, "doses")
  S <- attr(x, "S")
  off <- attr(x, "off")
  scal <- attr(x, "scal")
  ## calculate mean response at doses
  doseSeq <- seq(min(doses), max(doses), length=len)
  muMat <- getResp(altModels, doseSeq)

  if(is.null(x$sim))
    stop("Additional metrics only available if simulations were performed")
  cat("Calculating additional summary metrics, please wait.")
  ## calculate average mse of placebo-adjusted dose-response for ANOVA
  CM <- cbind(-1, diag(length(doses)-1))
  mseANOVA <- mean(diag(CM%*%S%*%t(CM)))
  ## calculate predictions
  predList <- getSimEst(x, "dose-response", doseSeq=doseSeq)
  
  out <- matrix(ncol = 5, nrow = ncol(muMat))
  colnames(out) <- c("Eff-vs-ANOVA", "cRMSE", "lengthTDCI", "P(no TD)", "lengthEDCI")
  rownames(out) <- colnames(muMat)
  if(!is.null(Delta)){
    direction <- match.arg(direction)
    tds <- getSimEst(x, "TD", Delta=Delta, direction=direction)
  }
  if(!is.null(p)){
    eds <- getSimEst(x, "ED", p=p)
  }
  for(i in 1:ncol(muMat)){
    out[i,1] <- mseANOVA/x$sim[i,1]^2
    ## calculate mse of estimating the plac-adj dose-response at fine grid
    ## first calculate placebo-adjusted predictions
    pred <- predList[[i]]
    pred <- (pred-pred[,1])[,-1]
    ## placebo-adjusted response
    mn <- (muMat[-1,i]-muMat[1,i])
    clmn <- colMeans(sweep(pred, 2, mn)^2)
    out[i,2] <- sqrt(mean(clmn))
    ## calculate length of CI for TD
    if(!is.null(Delta)){
      out[i,3] <- diff(quantile(tds[[i]], c(dLB, dUB), na.rm = TRUE))
      out[i,4] <- mean(is.na(tds[[i]]))
    } else {
      out[i,3] <- out[i,4] <- NA
    }
    ## calculate length of CI for ED
    if(!is.null(p)){
      out[i,5] <- diff(quantile(eds[[i]], c(dLB, dUB)))
    } else {
      out[i,5] <- NA
    }
  }
  cat("\r")
  cat(sprintf("Additional simulation metrics (nSim=%i)\n",
              attr(x$sim, "nSim")))
  print(signif(out, digits=digits))
}

## calculate the predictions for the fitted models
getSimEst <- function(x, type = c("dose-response", "ED", "TD"),
                      doseSeq, direction, p, Delta, placAdj = FALSE){
  modelSel <- attr(x$sim, "modelSel") 
  model <- attr(x, "model")
  coefs <- attr(x$sim, "coefs")
  off <- attr(x, "off")
  scal <- attr(x, "scal")
  nSim <- attr(x$sim, "nSim")
  altModels <- attr(x, "altModels")
  nAlt <- modCount(altModels, fullMod=TRUE)
  doses <- attr(x, "doses")
  maxD <- max(doses)
  type <- match.arg(type)
  if(type == "TD" & missing(direction))
    stop("need direction for TD calculation")

  out <- vector("list", nAlt)
  for(i in 1:nAlt){
    ind <- matrix(ncol = length(model), nrow = nSim)
    if(type == "dose-response"){
      resMat <- matrix(nrow = nSim, ncol = length(doseSeq))
      colnames(resMat) <- doseSeq
      rownames(resMat) <- 1:nSim
      for(j in 1:length(model)){
        ind[,j] <- modelSel[,i] == model[j]
        if(any(ind[,j])){
          cf <- do.call("rbind", (coefs[[i]])[ind[,j]])
          resMat[ind[,j]] <- predSamples(samples=cf, placAdjfullPars = TRUE,
                                         doseSeq=doseSeq,
                                         placAdj=placAdj, model=model[j],
                                         scal=scal, off=off, nodes = NULL)
        }
        out[[i]] <- resMat
      }
    }
    if(is.element(type, c("TD", "ED"))){
      resVec <- numeric(nSim)
      for(j in 1:length(model)){
        ind[,j] <- modelSel[,i] == model[j]
        if(any(ind[,j])){
          cf <- do.call("rbind", (coefs[[i]])[ind[,j]])
          if(type == "TD"){
            resVec[ind[,j]] <- apply(cf, 1, function(z){
              calcTD(model[j], z, Delta, "continuous", direction,
                     off=off, scal=scal)
            })
          }
          if(type == "ED"){
            resVec[ind[,j]] <- apply(cf, 1, function(z){
              calcED(model[j], z, p, maxD, "continuous", off=off, scal=scal)
            })
          }
        } 
      }
      out[[i]] <- resVec
    }
  }
  names(out) <- colnames(getResp(attr(x, "altModels"), doses=0)) ## horrible hack need to improve!
  out
}
  

plotDoseSims <- function(x, type = c("ED", "TD"), p, Delta, direction, xlab){
  altMods <- attr(x, "altModels")
  if(type == "ED"){
    out <- getSimEst(x, "ED", p=p)
    trueDoses <- ED(altMods, p=p, EDtype="continuous")
  } else {
    out <- getSimEst(x, "TD", Delta=Delta, direction=direction)
    trueDoses <- TD(altMods, Delta=Delta, TDtype="continuous",
                    direction=direction)
  }
  ## write plotting data frame
  nams <- names(out)
  group <- factor(rep(1:length(nams), each=length(out[[1]])), labels=nams)
  pdat <- data.frame(est = do.call("c", out),
                     group = group)
  ## determine limits for x-axis
  rngQ <- tapply(pdat$est, pdat$group, function(x){
    quantile(x, c(0.025, 0.975), na.rm=TRUE)
  })
  rngQ <- do.call("rbind", rngQ)
  rng <- c(min(rngQ[,1], na.rm = TRUE), max(rngQ[,2], na.rm = TRUE))
  delt <- diff(rng)*0.04
  ## truncate x-axis to 2*maxdose
  maxdose <- max(attr(x, "doses"))
  xlimU <- min(2*maxdose, max(rng[2], maxdose)+delt)
  xlimL <- max(-0.05*maxdose, min(0, rng[1])-delt)
  xlim <- c(xlimL, xlimU)
  parVal <- ifelse(type == "ED", paste("p=", p, sep=""), paste("Delta=", Delta, sep=""))
  maintxt <- paste("95%, 80%, 50% intervals and median of simulated ", type,
                   " estimates (", parVal, ")", sep = "")
  key <- list(text = list(maintxt, cex = 0.9))
  bwplot(~est|group, data=pdat, xlab = xlab, trueDoses=trueDoses,
         xlim = xlim,
         panel = function(...){
           z <- panel.number()
           panel.grid(v=-1, h=0, lty=2, col = "lightgrey")
           panel.abline(v=trueDoses[z], col = "red", lwd=2)
           panel.abline(v=c(0, max(attr(x, "doses"))), col = "grey", lwd=2)
           probs <- c(0.025, 0.1, 0.25, 0.5, 0.75, 0.9, 0.975)
           simDoseEst <- list(...)$x
           quants <- quantile(simDoseEst, probs, na.rm = TRUE)
           llines(c(quants[1], quants[7]), c(1,1), lwd=2, col=1)
           llines(c(quants[2], quants[6]), c(1,1), lwd=5, col=1)
           llines(c(quants[3], quants[5]), c(1,1), lwd=10, col=1)
           lpoints(quants[4], 1, cex=2, pch="|", col=1)
           if(type == "TD")
             ltext(xlim[2], 1.5, pos = 2, cex = 0.75,
                   labels = paste("% No TD:", mean(is.na(simDoseEst))*100, "%"))
         }, layout = c(1,length(out)), as.table = TRUE, key = key)
}

plotDRSims <- function(x, placAdj = FALSE, xlab, ylab){
  altMods <- attr(x, "altModels")
  rng <- range(attr(x, "doses"))
  doseSeq <- seq(rng[1], rng[2], length = 51)
  out <- getSimEst(x, type = "dose-response", doseSeq=doseSeq, placAdj = placAdj)
  trueMn <- getResp(altMods, doses=doseSeq)
  if(placAdj){
    trueMn <- trueMn-trueMn[1,]
  }
  nM <- length(out)
  resp <- vector("list", length=nM)
  for(i in 1:nM){
    qMat <-apply(out[[i]], 2, function(y){
      quantile(y, c(0.025, 0.25, 0.5, 0.75, 0.975))
    })
    resp[[i]] <- c(t(qMat))
  }
  
  resp <- do.call("c", resp)
  quant <- rep(rep(c(0.025, 0.25, 0.5, 0.75, 0.975), each = 51), nM)
  dose <- rep(doseSeq, nM*5)
  model <- factor(rep(1:nM, each = 5*51), labels = names(out))
  key <- list(text =
              list("Pointwise 95%, 50% intervals and median of simulated dose-response estimates", cex = 0.9))

  xyplot(resp~dose|model, groups = quant, xlab=xlab, ylab = ylab,
         panel = function(...){
           ## plot grid
           panel.grid(v=-1, h=-1, col = "lightgrey", lty=2)
           ## plot estimates
           panel.dat <- list(...)
           ind <- panel.dat$subscripts
           LB95.x <- panel.dat$x[panel.dat$groups[ind] == 0.025]
           LB95 <- panel.dat$y[panel.dat$groups[ind] == 0.025]
           UB95.x <- panel.dat$x[panel.dat$groups[ind] == 0.975]
           UB95 <- panel.dat$y[panel.dat$groups[ind] == 0.975]
           lpolygon(c(LB95.x, rev(UB95.x)), c(LB95, rev(UB95)),
                    col = "lightgrey", border = "lightgrey")
           LB50.x <- panel.dat$x[panel.dat$groups[ind] == 0.25]
           LB50 <- panel.dat$y[panel.dat$groups[ind] == 0.25]
           UB50.x <- panel.dat$x[panel.dat$groups[ind] == 0.75]
           UB50 <- panel.dat$y[panel.dat$groups[ind] == 0.75]
           lpolygon(c(LB50.x, rev(UB50.x)), c(LB50, rev(UB50)),
             col = "darkgrey", border = "darkgrey")
           MED.x <- panel.dat$x[panel.dat$groups[ind] == 0.5]
           MED <- panel.dat$y[panel.dat$groups[ind] == 0.5]
           llines(MED.x, MED, col = 1,lwd = 1.5)
           ## plot true curve
           z <- panel.number()
           llines(doseSeq, trueMn[,z], col=2, lwd=1.5)
         }, as.table = TRUE, key=key)
}


plot.planMod <- function(x, type = c("dose-response", "ED", "TD"),
                         p, Delta, direction, placAdj = FALSE,
                         xlab = "Dose", ylab = "", ...){
  type <- match.arg(type)
  if(type == "dose-response"){
    plotDRSims(x, placAdj = placAdj, xlab=xlab, ylab = ylab)
  } else {
    plotDoseSims(x, type=type, p=p, Delta=Delta,
                 direction = direction, xlab = xlab)
  }
}
