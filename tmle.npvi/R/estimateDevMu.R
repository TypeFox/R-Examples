estimateDevMu <- function(muW, obs, eic1, flavor=c("learning", "superLearning"), learnDevMu,
                          light=TRUE, SuperLearner.=NULL, ..., verbose=FALSE) {
  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ## Validate arguments
  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ## Argument 'mu':
  muW <- Arguments$getNumerics(muW);
  
  ## Argument 'obs':
  obs <- validateArgumentObs(obs, allowIntegers=TRUE);
  
  ## Argument 'eic1'
  eic1 <- Arguments$getNumerics(eic1);
  
  ## Argument 'flavor':
  flavor <- match.arg(flavor);
  learnDevMode <- switch(flavor,
                      learning="function",
                      superLearning="character");

  ## Argument 'learnDevMu'
  mode <- mode(learnDevMu);
  if (mode != learnDevMode) {
    throw("Argument 'learnDevMu' should be of mode '", learnDevMode, "', not '", mode, "' for flavor: ", flavor);
  }

  ## Argument 'SuperLearner.'
  if (flavor=="superLearning") {
    if (is.null(SuperLearner.) || mode(SuperLearner.)!="function") {
      throw("Argument 'SuperLearner.' should be a function")
    }
  }

  
  ## Argument 'verbose'
  verbose <- Arguments$getVerbose(verbose);

  if (flavor=="learning") {
    devMu <- learnDevMu(obs, eic1, muW, light=light, ...);
  } else if (flavor=="superLearning") {
    logSL <- as.logical(less(verbose, 10));  ## decrease verbosity in superLearner
    SL.library.devMu <- learnDevMu
    obsD <- as.data.frame(obs)
    ZdevMu <- (obsD[, "X"] - muW) * eic1;

    fitDevMu <- SuperLearner.(Y=ZdevMu, X=extractW(obsD), ## obsD[, "W", drop=FALSE]
                              SL.library=SL.library.devMu, verbose=logSL,
                              family=gaussian(), ...);
    devMu <- function(W) {
      Wd <- as.data.frame(W)
      predict(fitDevMu, newdata=Wd)$pred;
    }
  }
  verbose && cat(verbose, "devMu(W):");
  verbose && print(verbose, summary(devMu(extractW(obs))));
  
  devMu
}


############################################################################
## HISTORY:
## 2014-02-07
## o Created.
############################################################################

