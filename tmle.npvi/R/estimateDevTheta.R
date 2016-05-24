estimateDevTheta <- function(thetaXW, obs, flavor=c("learning", "superLearning"), learnDevTheta,
                             light=TRUE, SuperLearner.=NULL, ..., verbose=FALSE) {
  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ## Validate arguments
  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ## Argument 'thetaXW':
  thetaXW <- Arguments$getNumerics(thetaXW);
    
  ## Argument 'obs':
  obs <- validateArgumentObs(obs, allowIntegers=TRUE);

  ## Argument 'flavor':
  flavor <- match.arg(flavor);
  learnDevMode <- switch(flavor,
                      learning="function",
                      superLearning="character");

  ## Argument 'learnDevTheta'
  mode <- mode(learnDevTheta);
  if (mode != learnDevMode) {
    throw("Argument 'learnDevTheta' should be of mode '", learnDevMode, "', not '", mode, "' for flavor: ", flavor);
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
    devTheta <- learnDevTheta(obs, thetaXW, light=light, verbose=verbose);
  } else if (flavor=="superLearning") {
    logSL <- as.logical(less(verbose, 10));  ## decrease verbosity in superLearner
    obsD <- as.data.frame(obs)
    ZdevTheta <- (obsD[, "Y"]-thetaXW)^2;
    SL.library.devTheta <- learnDevTheta;

    fitDevTheta <- SuperLearner.(Y=ZdevTheta, X=extractXW(obsD),  ## obsD[, c("X", "W")]
                                 SL.library=SL.library.devTheta, verbose=logSL,
                                 family=gaussian(), ...);
    
    devTheta <- function(XW) {
      XWd <- as.data.frame(XW)
      predict(fitDevTheta, newdata=XWd)$pred
    }
  }
  verbose && cat(verbose, "devTheta(XW):");
  verbose && print(verbose, summary(devTheta(extractXW(obs))));

  devTheta
}


############################################################################
## HISTORY:
## 2014-02-07
## o Created.
############################################################################

