estimateDevG <- function(gW, obs, eic1, flavor=c("learning", "superLearning"), learnDevG,
                         light=TRUE, SuperLearner.=NULL, ..., verbose=FALSE) {
  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ## Validate arguments
  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ## Argument 'gW':
  gW <- Arguments$getNumerics(gW);
  
  ## Argument 'obs':
  obs <- validateArgumentObs(obs, allowIntegers=FALSE);

  ## Argument 'eic1'
  eic1 <- Arguments$getNumerics(eic1);
  
  ## Argument 'flavor':
  flavor <- match.arg(flavor);
  learnDevMode <- switch(flavor,
                      learning="function",
                      superLearning="character");

  ## Argument 'learnDevG'
  mode <- mode(learnDevG);
  if (mode != learnDevMode) {
    throw("Argument 'learnDevG' should be of mode '", learnDevMode, "', not '", mode, "' for flavor: ", flavor);
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
    devG <- learnDevG(obs, eic1, gW, light=light, ...);
  } else if (flavor=="superLearning") {
    logSL <- as.logical(less(verbose, 10));  ## decrease verbosity in SuperLearner
    obsD <- as.data.frame(obs)
    ZdevG <- eic1 * ( (obsD[, "X"]==0) - gW );
    SL.library.devG <- learnDevG;

    fitDevG <- SuperLearner.(Y=ZdevG, X=extractW(obsD), ## obsD[, "W", drop=FALSE]
                             SL.library=SL.library.devG, verbose=logSL,
                             family=gaussian(), ...)
    devG <- function(W) {
      Wd <- as.data.frame(W)
      predict(fitDevG, newdata=Wd)$pred
    }
  }
  verbose && cat(verbose, "devG(W):");
  verbose && print(verbose, summary(devG(extractW(obs))));

  devG
}


############################################################################
## HISTORY:
## 2014-02-07
## o Created.
############################################################################

