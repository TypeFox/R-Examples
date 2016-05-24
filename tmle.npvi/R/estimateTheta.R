estimateTheta <- function(obs, flavor=c("learning", "superLearning"), learnTheta,
                          light=TRUE, SuperLearner.=NULL, ..., verbose=FALSE) {
  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ## Validate arguments
  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ## Argument 'obs':
  obs <- validateArgumentObs(obs);

  ## Argument 'flavor':
  flavor <- match.arg(flavor);
  learnMode <- switch(flavor,
                      learning="function",
                      superLearning="character");

  ## Argument 'learnTheta'
  mode <- mode(learnTheta);
  if (mode != learnMode) {
    throw("Argument 'learnTheta' should be of mode '", learnMode, "', not '", mode, "' for flavor: ", flavor);
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
    theta <- learnTheta(obs, light=light, ...);
  } else if (flavor=="superLearning") {
    logSL <- as.logical(less(verbose, 10));  ## decrease verbosity in SuperLearner
    SL.library.theta <- learnTheta;
    obsD <- as.data.frame(obs)

    fitTheta <- SuperLearner.(Y=obsD[, "Y"], X=extractXW(obsD), ## obsD[, c("X", "W")]
                              SL.library=SL.library.theta, verbose=logSL,
                              family=gaussian(), ...)
    theta <- function(XW) {
      XWd <- as.data.frame(XW)
      predict(fitTheta, newdata=XWd)$pred
    }
  }
  verbose && cat(verbose, "theta(X,W):");
  verbose && print(verbose, summary(theta(extractXW(obs)))); ## obs[, c("X", "W")]

  theta
}


############################################################################
## HISTORY:
## 2014-02-07
## o Created.
############################################################################

