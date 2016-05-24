setMethodS3("getMumin", "NPVI", function(this, ...) {
  this$.mumin;
})

setMethodS3("getMumax", "NPVI", function(this, ...) {
  this$.mumax;
})

setMethodS3("getMu", "NPVI", function(this, tabulate, ...) {
  if (missing(tabulate)) {
    tabulate <- getTabulate(this);
  }
  if (!tabulate) {
    this$.mu;
  } else {
    this$.mutab;
  } 
})

setMethodS3("getMuAux", "NPVI", function(this, tabulate, ...) {
  if (missing(tabulate)) {
    tabulate <- getTabulate(this);
  }
  if (!tabulate) {
    this$.muAux;
  } else {
    this$.muAuxtab;
  } 
})


setMethodS3("setMu", "NPVI", function(this, mu, ...) {
  ## Argument 'mu':
  if ((!is.null(mu))  && (mode(mu)!="function")) {
    throw("Argument \var{mu} should be of mode 'function', not ", mode(mu));
  }
  
  mumin <- getMumin(this)
  mumax <- getMumax(this)
  thresholdedMu <- function(W) {
    threshold(mu(W), min=mumin, max=mumax)
  }
  this$.mu <- thresholdedMu ;
})

setMethodS3("setMuAux", "NPVI", function(this, muAux, ...) {
  ## Argument 'muAux':
  if ((!is.null(muAux))  && (mode(muAux)!="function")) {
    throw("Argument \var{muAux} should be of mode 'function', not ", mode(muAux));
  }
  
  mumin <- getMumin(this)
  mumax <- getMumax(this)
  thresholdedMuAux <- function(W) {
    threshold(muAux(W), min=mumin, max=mumax)
  }
  this$.muAux <- thresholdedMuAux ;
})


setMethodS3("setMuTab", "NPVI", function(this, mu, ...) {
  ## Argument 'mu':
  if ((!is.null(mu))  && (mode(mu)!="function")) {
    throw("Argument \var{mu} should be of mode 'function', not ", mode(mu));
  }

  mumin <- getMumin(this)
  mumax <- getMumax(this)
  thresholdedMu <- function(W) {
    threshold(mu(W), min=mumin, max=mumax)
  }

  this$.mutab <- thresholdedMu
})

setMethodS3("setMuAuxTab", "NPVI", function(this, muAux, ...) {
  ## Argument 'muAux':
  if ((!is.null(muAux))  && (mode(muAux)!="function")) {
    throw("Argument \var{muAux} should be of mode 'function', not ", mode(muAux));
  }

  mumin <- getMumin(this)
  mumax <- getMumax(this)
  thresholdedMuAux <- function(W) {
    threshold(muAux(W), min=mumin, max=mumax)
  }

  this$.muAuxtab <- thresholdedMuAux
})


setMethodS3("initializeMu", "NPVI", function(this, muAux, g, ...) {
  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ## Validate arguments
  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ## Argument 'muAux':
  if (mode(muAux) != "function") {
    throw("Argument 'muAux' should be a function, not a ", mode(muAux));
  }
  ## Argument 'g':
  if (mode(g) != "function") {
    throw("Argument 'g' should be a function, not a ", mode(g));
  }

  setMuAux(this, muAux)

  ## mu
  mu <- function(W) {
    muAux(W)*(1-g(W))
  }
  setMu(this, mu)

  ## tabulated version of 'muAux'
  fW <- getFW(this);
  obs <- getObs(this);
  MUAUXTAB <- muAux(fW(obs)); ## a *vector*, not a function
  muAuxtab <- function(ii) {
    MUAUXTAB[ii];
  }
  setMuAuxTab(this, muAuxtab)
  
  ## tabulated version of 'mu'
  MUTAB <- mu(fW(obs)); ## a *vector*, not a function
  mutab <- function(ii) {
    MUTAB[ii];
  }
  setMuTab(this, mutab)
})

setMethodS3("updateMu", "NPVI", function(this, dev, exact=TRUE, effICW, ...) {
  updateMuNonTab(this, dev, exact=exact, effICW, ...)
  updateMuTab(this, dev, exact=exact, effICW, ...)
})

setMethodS3("updateMuNonTab", "NPVI", function(this, dev, exact=TRUE, effICW, ...) {
  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ## Validate arguments
  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ## Argument 'dev':
  if (mode(dev) != "function") {
    throw("Argument 'dev' should be a function, not a ", mode(dev));
  }
  ## Argument 'exact':
  exact <- Arguments$getLogical(exact);
  ## Argument 'effICW':
  if (exact) {
    if (missing(effICW)) {
      throw("Argument 'effICW' is required when 'exact' is TRUE");
    }
    if (mode(effICW) != "function") {
      throw("Argument 'dev' should be a function, not a ", mode(dev));
    }
  }

  fW <- getFW(this)
  eps <- getEpsilon(this)

  mu <- getMu(this, tabulate=FALSE);
  g <- getG(this, tabulate=FALSE);
  
  if (!exact) { ## if do not use exact expression
    mu1 <- function(W) {
      mu(W) + eps * dev(W);
    }
  } else { ## if use exact expression
    mu1 <- function(W) {
      muW <- mu(W);
      theEffICW <- effICW(W)
      numerator <- muW + eps * (dev(W) + muW*theEffICW);
      denominator <- 1 + eps*theEffICW;
      numerator/denominator;
    }
  }
  muAux1 <- function(W) {
    mu1(W)/(1-g(W))
  }
  setMuAux(this, muAux1);
  setMu(this, mu1);
})

setMethodS3("updateMuTab", "NPVI", function(this, dev, exact=TRUE, effICW, ...) {
  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ## Validate arguments
  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ## Argument 'dev':
  if (mode(dev) != "function") {
    throw("Argument 'dev' should be a function, not a ", mode(dev));
  }
  ## Argument 'exact':
  exact <- Arguments$getLogical(exact);
  
  ## Argument 'effICW':
  if (exact) {
    if (missing(effICW)) {
      throw("Argument 'effICW' is required when 'exact' is TRUE");
    }
    if (mode(effICW) != "function") {
      throw("Argument 'dev' should be a function, not a ", mode(dev));
    }
  }

  fW <- getFW(this)
  eps <- getEpsilon(this)

  mu <- getMu(this, tabulate=TRUE)
  g <- getG(this, tabulate=TRUE)
  obs <- getObs(this, tabulate=TRUE);
  muW <- mu(obs[, "W"])
  gW <- g(obs[, "W"])
  rm(mu, g, obs)

  obs <- getObs(this)
  devW <- dev(fW(obs))
  W <- obs[, "W"]
  rm(obs)
  
  if (!exact) { ## do not use exact expression
    mu1W <- muW + eps * devW;
  } else { ## use exact expression
    theEffICW <- effICW(W)
    ## the above should use the tabulated or real versions of mu and theta0
    ## depending on tabulate, because effICW works on true values or
    ## indices depending on 'tabulate' (see how it is created in 'NPVI.update')
    numerator <- muW + eps * (devW + muW*theEffICW);
    denominator <- 1 + eps*theEffICW;
    mu1W <- numerator/denominator;
  }

  muAux1W <- mu1W/(1-gW)
  muAux1tab <- function(ii) {
    muAux1W[ii]
  }
  setMuAuxTab(this, muAux1tab)

  mu1tab <- function(ii) {
    mu1W[ii]
  }
  setMuTab(this, mu1tab)
})

############################################################################
## HISTORY:
## 2014-02-07
## o Created.
############################################################################

