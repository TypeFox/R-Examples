estimateMuAux <- function(obs, flavor=c("learning", "superLearning"), learnMuAux,
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

  ## Argument 'learnMuAux'
  mode <- mode(learnMuAux);
  if (mode != learnMode) {
    throw("Argument 'learnMuAux' should be of mode '", learnMode, "', not '", mode, "' for flavor: ", flavor);
  }

  ## Argument 'SuperLearner.'
  if (flavor=="superLearning") {
    if (is.null(SuperLearner.) || mode(SuperLearner.)!="function") {
      throw("Argument 'SuperLearner.' should be a function")
    }
  }

  
  ## Argument 'verbose'
  verbose <- Arguments$getVerbose(verbose);

  idx <- which(obs[, "X"] != 0);
  
  if (flavor=="learning") {
    muAux <- learnMuAux(obs[idx, ], light=light, ...);
  } else if (flavor=="superLearning") {
    logSL <- as.logical(less(verbose, 10));  ## decrease verbosity in SuperLearner
    SL.library.muAux <- learnMuAux;
    obsD <- as.data.frame(obs)
    ## WW <- obsD[idx, "W", drop=FALSE]
    WW <- extractW(obsD[idx, ])

    fitMuAux <- SuperLearner.(Y=obsD[idx, "X"], X=WW,
                              SL.library=SL.library.muAux, verbose=logSL,
                              family=gaussian(), ...);
    verbose && print(verbose, fitMuAux);
    muAux <- function(W) {
      Wd <- as.data.frame(W)
      predict(fitMuAux, newdata=Wd)$pred;
    }
  }
  verbose && cat(verbose, "mu'(W):");
  verbose && print(verbose, summary(muAux(extractW(obs))));
  
  muAux
}

############################################################################
## HISTORY:
## 2014-02-07
## o Created.
############################################################################

