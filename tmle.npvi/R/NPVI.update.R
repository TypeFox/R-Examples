setMethodS3("update", "NPVI", function(object,
                                       flavor=c("learning", "superLearning"),
                                       cvControl=NULL,
                                       learnDevG=NULL,
                                       learnDevMu=NULL,
                                       learnDevTheta=NULL,
                                       learnCondExpX2givenW=NULL,
                                       learnCondExpXYgivenW=NULL,
                                       bound=1e-1, B=1e4,
                                       light=TRUE, 
                                       cleverCovTheta=TRUE,
                                       exact=TRUE, trueGMu=NULL,
                                       SuperLearner.=SuperLearner.,
                                       ..., verbose=FALSE) {

  this <- object; ## to please R CMD CHECK 
  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ## Validate arguments
  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ## Argument 'flavor'
  flavor <- match.arg(flavor);
  learnDevMode <- switch(flavor,
                         learning="function",
                         superLearning="character");
  
  ## Argument 'learnDevG'
  mode <- mode(learnDevG);
  if (mode != learnDevMode) {
    throw("Argument 'learnDevG' should be of mode '", learnDevMode, "', not '", mode, "' for flavor: ", flavor);
  }

  ## Argument 'learnDevMu'
  mode <- mode(learnDevMu);
  if (mode != learnDevMode) {
    throw("Argument 'learnDevMu' should be of mode '", learnDevMode, "', not '", mode, "' for flavor: ", flavor);
  }

  ## Argument 'learnDevTheta'
  mode <- mode(learnDevTheta);
  if ((mode != learnDevMode) & (!cleverCovTheta)) {
    throw("Argument 'learnDevTheta' should be of mode '", learnDevMode, "', not '", mode, "' for flavor: ", flavor);
  }

  ## Argument 'bound':
  bound <- Arguments$getNumeric(bound);
  if (bound<=0) {
    throw("Argument 'bound' must be positive!\n")
  }
  ## Argument 'B':
  B <- Arguments$getInteger(B);

  ## Argument 'light'
  light <- Arguments$getLogical(light);

  ## Argument 'cleverCovTheta'
  cleverCovTheta <- Arguments$getLogical(cleverCovTheta);

  ## Argument 'exact'
  exact <- Arguments$getLogical(exact);

  ## Argument 'trueGMu'
  useTrueGMu <- (!is.null(trueGMu))
  
  ## Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  verbose <- less(verbose, 10);

  ## Incrementing the current value of 'step'
  this$.step <- getStep(this)+1;

  ## Argument 'SuperLearner.'
  if (flavor=="superLearning") {
    if (is.null(SuperLearner.) || mode(SuperLearner.)!="function") {
      throw("Argument 'SuperLearner.' should be a function")
    }
  }

  
  verbose && cat(verbose, "Iteration ", getStep(this), "\n");


  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ## Retrieve elements of 'this'
  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  obs <- getObs(this);
  obsT <- getObs(this, tabulate=FALSE);
  Xq <- getXq(this);
  Yq <- getYq(this);
  family <- getFamily(this);
  tabulate <- getTabulate(this)
  g <- getG(this);
  mu <- getMu(this);
  muAux <- getMuAux(this);
  theta <- getTheta(this);
  theta0 <- getTheta0(this);
  sigma2 <- getSigma2(this);
  psi <- getPsi(this);
  fW <- getFW(this)
  fX <- getFX(this)
  
  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ## Update divergence (part 1/2)
  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (cleverCovTheta) {
    div <- NA  ## 'div' cannot be calculated without further (otherwise unnecessary) assumptions
  } else {
    weightsW <- getWeightsW(this)
    fY <- getFY(this)
    obsB <- simulateData(B, obs[, "W"], obsT[, "X"], Xq, g, mu, muAux, sigma2,
                         theta=theta, Y=Yq, 
                         weightsW=weightsW,
                         family=family)
    ## taken from 'updateEfficientInfluenceCurve'
    thetaXW <- theta(obsB[, c("X", "W")]);
    theta0W <- theta0(obsB[, "W", drop=FALSE]);
    muW <- mu(obsB[, "W", drop=FALSE])
    gW <- g(obsB[, "W", drop = FALSE])
    
    X <- fX(obsB)
    Y <- fY(obsB)
              
    D1 <- X * (thetaXW - theta0W - X * psi);
    D2 <- (Y - thetaXW) * (X - muW/gW*(X==0));
    verbose && summary(verbose, D1);
    verbose && summary(verbose, D2);
    verbose && print(verbose, sigma2);
    
    eic1 <- D1 / sigma2;
    eic2 <- D2 / sigma2;
    eic <- eic1 + eic2;
    verbose && summary(verbose, eic);

    partialDiv <- mean(abs(eic))
    rm(eic1, eic2, eic);
  }
  
  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ## Efficient influence curve
  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  eic <- getEfficientInfluenceCurve(this);
  eic1 <- eic[, "eic1"];
  rm(eic);

  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ## Update the estimation of relevant components of the distribution
  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Updating the estimation of relevant components of the distribution");

  if (cleverCovTheta) {
    updateEpsilonTheta(this);
    ## Update 'theta', then update EIC and 'epsilon' accordinlgy
    updateTheta(this, NULL, cleverCovTheta=cleverCovTheta, exact=exact);
    updatePsi(this, B, verbose=verbose);
    updateEfficientInfluenceCurve(this, obs);
    
    theta <- getTheta(this)
    theta0 <- getTheta0(this)
    psi <- getPsi(this)
    eic <- getEfficientInfluenceCurve(this)
    eic1 <- eic[, "eic1"]
    rm(eic)

    ## Update history
    setDivergence(this, div)
    updateHistory(this);
    if (!useTrueGMu) {
      ## browser(skipCalls=2)
      updateEpsilon(this, cleverCovTheta=TRUE, bound=bound);  ## preparing update of 'mu' and 'g'
    }
  } else {  ## update 'theta' without using a clever covariate, and don't update EIC and 'epsilon'
    updateEpsilon(this, cleverCovTheta=FALSE, bound=bound);
    div <- abs(getEpsilon(this))*partialDiv
    thetaXW <- theta(obs[, c("X", "W")])
    devTheta <- estimateDevTheta(thetaXW, obsT, flavor=flavor, learnDevTheta=learnDevTheta, light=light,
                                 SuperLearner.=SuperLearner., ..., verbose=verbose);
    updateTheta(this, devTheta, cleverCovTheta=cleverCovTheta, exact=exact);
  }

  effICW <- NULL;

  if (!useTrueGMu) {
    ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ## Update estimation of 'g' and 'mu'
    ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    if (flavor=="learning") {
      condExpX2givenW <- learnCondExpX2givenW(obsT, light=light); ## a 'true' function
      condExpXYgivenW <- learnCondExpXYgivenW(obsT, light=light); ## a 'true' function
    } else if (flavor=="superLearning") {
      logSL <- as.logical(less(verbose, 10));  ## decrease verbosity in SuperLearner
      SL.library.condExpX2givenW <- learnCondExpX2givenW; 
      SL.library.condExpXYgivenW <- learnCondExpXYgivenW; 
      obsD <- as.data.frame(obsT)
      fitCondExpX2givenW <- SuperLearner.(Y=obsD[, "X"]^2, X=extractW(obsD), ## obsD[, "W", drop=FALSE],
                                          SL.library=SL.library.condExpX2givenW, verbose=logSL,
                                          family=gaussian(), ...)
      fitCondExpXYgivenW <- SuperLearner.(Y=obsD[, "Y"]*obsD[, "X"], X=extractW(obsD), ## obsD[, "W", drop=FALSE],
                                          SL.library=SL.library.condExpXYgivenW, verbose=logSL,
                                          family=gaussian(), ...)
      condExpX2givenW <- function(W) {
        Wd <- as.data.frame(W)
        predict(fitCondExpX2givenW, newdata=Wd)$pred
      }
      condExpXYgivenW <- function(W) {
        Wd <- as.data.frame(W)
        predict(fitCondExpXYgivenW, newdata=Wd)$pred
      }
      verbose && cat(verbose, "E(X^2|W):");
      verbose && print(verbose, summary(condExpX2givenW(extractW(obsD))));
      ##
      verbose && cat(verbose, "E(XY|W):");
      verbose && print(verbose, summary(condExpXYgivenW(extractW(obsD))));
    }
    
    ## pasted from 'estimateEpsilon'
    theMin <- min(eic1)
    theMax <- max(eic1)
    if (theMin > 0) {
      interval <- c(-1/(1.001*theMax), 1e3)
    }
    if (theMax < 0) {
      interval <- c(-1e3, -1/(1.001*theMin))
    }
    if (theMin<=0 & theMax>=0) {
      interval <- c(-1/(1.001*theMax), -1/(1.001*theMin))
    }
   
    effICW <- function(W) {
      realW <- fW(cbind(W=W, X=NA, Y=NA))
      out <- (condExpXYgivenW(realW)-theta0(W)*mu(W)-psi*condExpX2givenW(realW))/sigma2;
      threshold(out, theMin, theMax)
    }
    
    
    ## Update 'g' *before* 'mu' as the updated 'mu' depends on the updated 'g', see *inside* 'updateMu'.
    gW <- g(extractW(obs))
    devG <- estimateDevG(gW, obsT, eic1, flavor=flavor, learnDevG=learnDevG, light=light,
                         SuperLearner.=SuperLearner.,
                         ..., verbose=verbose);
    updateG(this, devG, exact=exact, effICW=effICW);


    muW <- mu(extractW(obs))
    devMu <- estimateDevMu(muW, obsT, eic1, flavor=flavor, learnDevMu=learnDevMu, light=light,
                           SuperLearner.=SuperLearner.,
                           ..., verbose=verbose);
    updateMu(this, devMu, exact=exact, effICW=effICW);


    ## Update 'sigma2'
    X <- fX(obs)
    devSigma2 <- mean(eic1 * X^2);
    updateSigma2(this, devSigma2);  
    verbose && exit(verbose);
    updateWeightsW(this, effICW)    
  }

  
  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ## Update estimation of 'psi'
  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Updating the estimation of 'psi'");
  updatePsi(this, B, verbose=verbose);
  psi <- getPsi(this);
  verbose && cat(verbose, "psi: ", round(psi, 3));

  updateEfficientInfluenceCurve(this);  
  
  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ## Update history
  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  setDivergence(this, div)
  updateHistory(this);

  verbose && exit(verbose);
})

############################################################################
## HISTORY:
## 2014-02-07
## o Created.
############################################################################

