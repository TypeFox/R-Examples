setMethodS3("getThetamin", "NPVI", function(this, ...) {
  this$.thetamin;
})

setMethodS3("getThetamax", "NPVI", function(this, ...) {
  this$.thetamax;
})

setMethodS3("getTheta", "NPVI", function(this, tabulate, ...) {
  if (missing(tabulate)) {
    tabulate <- getTabulate(this);
  }
  if (!tabulate) {
    this$.theta;
  } else {
    this$.thetatab;
  } 
})

setMethodS3("setTheta", "NPVI", function(this, theta, ...) {
  ## Argument 'theta':
  if ((!is.null(theta))  && (mode(theta)!="function")) {
    throw("Argument \var{theta} should be of mode 'function', not ", mode(theta));
  }
  
  thetamin <- getThetamin(this)
  thetamax <- getThetamax(this)

  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ## set 'theta'
  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  thresholdedTheta <- function(XW) {
    threshold(theta(XW), min=thetamin, max=thetamax)
  }
  this$.theta <- thresholdedTheta;

  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ## set 'theta0' accordingly
  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  thresholdedTheta0 <- function(W) {
    XW <- cbind(X=0, W=W);
    threshold(theta(XW), min=thetamin, max=thetamax);
  }
  this$.theta0 <- thresholdedTheta0;
})

setMethodS3("setThetaTab", "NPVI", function(this, theta, theta0, ...) {
  ## Argument 'theta':
  if ((!is.null(theta))  && (mode(theta)!="function")) {
    throw("Argument \var{theta} should be of mode 'function', not ", mode(theta));
  }

  ## Argument 'theta0':
  if ((!is.null(theta0))  && (mode(theta0)!="function")) {
    throw("Argument \var{theta0} should be of mode 'function', not ", mode(theta0));
  }
  
  thetamin <- getThetamin(this)
  thetamax <- getThetamax(this)

  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ## set 'theta'
  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  thresholdedTheta <- function(XW) {
    threshold(theta(XW), min=thetamin, max=thetamax)
  }
  this$.thetatab <- thresholdedTheta

  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ## set 'theta0'
  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  thresholdedTheta0 <- function(W) {
    threshold(theta0(W), min=thetamin, max=thetamax)
  }
  this$.theta0tab <- thresholdedTheta0
})

setMethodS3("getTheta0", "NPVI", function(this, tabulate, ...) {
  if (missing(tabulate)) {
    tabulate <- getTabulate(this);
  }
  if (!tabulate) {
    this$.theta0;
  } else {
    this$.theta0tab;
  }
})

setMethodS3("initializeTheta", "NPVI", function(this, theta, ...) {
  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ## Validate arguments
  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ## Argument 'theta':
  if (mode(theta) != "function") {
    throw("Argument 'theta' should be a function, not a ", mode(theta));
  }

  ## theta
  setTheta(this, theta)

  ## tabulated versions of 'theta' and 'theta0'
  fW <- getFW(this);
  fX <- getFX(this);
  obs <- getObs(this);
  nr <- nrow(obs)
  Xq <- getXq(this);
  Xq.value <- Arguments$getNumerics(Xq$value)
  Xq.index <- Arguments$getIntegers(Xq$index)
  
  eg.value <- expand.grid(X=Xq.value, W=obs[, "W"])
  eg.index <- expand.grid(X=Xq.index, W=1:nr)
  ## non-diagonal entries
  OBSTAB.nondiag <- cbind(X=eg.value[, "X"], fW(eg.value))
  THETATAB.nondiag <- theta(OBSTAB.nondiag) ## length(Xq.value)*nrow(obs) vector
  ## diagonal entries
  OBSTAB.diag <- cbind(X=fX(obs), W=fW(obs))
  THETATAB.diag <- theta(OBSTAB.diag) ## a nrow(obs) vector
  ## making the sparse matrix
  THETATAB <- sparseMatrix(i=eg.index[, "X"], j=eg.index[, "W"],
                           x=as.vector(THETATAB.nondiag),
                           dims=c(nr, nr))
  THETATAB[cbind(1:nr, 1:nr)] <- as.vector(THETATAB.diag)
  THETA0TAB <- theta(cbind(X=0, W=fW(obs)))  ## a vector, not a matrix!
  thetatab <- function(xiwj) {
    stopifnot(is.matrix(xiwj) && ncol(xiwj)==2 && is.integer(xiwj))
    THETATAB[xiwj]
  }
  theta0tab <- function(wi) {
    ## stopifnot(is.matrix(wi) && ncol(wi)==1 && is.integer(wi))
    THETA0TAB[wi]
  }
  setThetaTab(this, thetatab, theta0tab);
})

setMethodS3("updateTheta", "NPVI", function(this, dev, cleverCovTheta, exact=TRUE, ...) {
  updateThetaNonTab(this, dev, cleverCovTheta, exact=exact, ...)
  updateThetaTab(this, dev, cleverCovTheta, exact=exact, ...)
})

setMethodS3("updateThetaNonTab", "NPVI", function(this, dev, cleverCovTheta, exact=TRUE, ...) {
  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ## Validate arguments
  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ## Argument 'dev':
  if ((mode(dev) != "function") & (!cleverCovTheta)) {
    throw("Argument 'dev' should be a function, not a ", mode(dev));
  }

  ## Argument 'cleverCovTheta':
  cleverCovTheta <- Arguments$getLogical(cleverCovTheta);

  ## Argument 'exact':
  exact <- Arguments$getLogical(exact);

  fW <- getFW(this)
  fX <- getFX(this)
  devThetaAux <- dev

  theta <- getTheta(this, tabulate=FALSE);
  
  if (!cleverCovTheta) {
    eps <- getEpsilon(this);
    g <- getG(this, tabulate=FALSE);
    mu <- getMu(this, tabulate=FALSE);
    sigma2 <- getSigma2(this);
    psi <- getPsi(this);
    if (!exact) { ## do not use clever covariate nor exact expression
      dev <- function(XW) {
        X <- XW[, 1, drop=TRUE];
        W <- XW[, 2, drop=FALSE];
        devThetaAux(XW) * (X - mu(W)/g(W)*(X==0))/sigma2;
      }
      theta1 <- function(XW) {
        theta(XW) + eps * dev(XW);
      }
    } else { ## do not use clever covariate, but exact expression
      theta1 <- function(XW) {
        X <- XW[, 1, drop=TRUE];
        W <- XW[, 2, drop=FALSE];
        TXW <- theta(XW);
        term1 <- devThetaAux(XW) * (X - mu(W)/g(W)*(X==0))/sigma2;
        term2 <- (X*(TXW-theta(cbind(X=0, W=W))-X*psi))/sigma2;
        numerator <- TXW + eps*(term1+term2*TXW);
        denominator <- 1 + eps*term2;
        numerator/denominator;
      }
    }
  } else { ## use clever covariate (hence, exact expression)
    if (!exact) {
      throw("Can't use a clever covariate for updating 'theta' if 'exact==FALSE' !");
    }
    eps <- getEpsilonTheta(this);
    dev <- getHTheta(this, tabulate=FALSE);
    theta1 <- function(XW) {
      theta(XW) + eps * dev(XW);
    }
  }
  setTheta(this, theta1)
})

setMethodS3("updateThetaTab", "NPVI", function(this, dev, cleverCovTheta, exact=TRUE, ...) {
  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ## Validate arguments
  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ## Argument 'dev':
  if ((mode(dev) != "function") & (!cleverCovTheta)) {
    throw("Argument 'dev' should be a function, not a ", mode(dev));
  }

  
  ## Argument 'cleverCovTheta':
  cleverCovTheta <- Arguments$getLogical(cleverCovTheta);

  ## Argument 'exact':
  exact <- Arguments$getLogical(exact);

  Xq <- getXq(this);
  Xq.value <- Arguments$getNumerics(Xq$value)
  Xq.index <- Arguments$getIntegers(Xq$index)
  
  fW <- getFW(this)
  fX <- getFX(this)
  devThetaAux <- dev

  theta <- getTheta(this, tabulate=TRUE);
  theta0 <- getTheta0(this, tabulate=TRUE)

  obs <- getObs(this, tabulate=TRUE);
  nr <- nrow(obs)

  ## very old:
  ## rm(obs)
  ## XW <- as.matrix(expand.grid(X=1:nr, W=1:nr))
  ## thetaXW <- matrix(theta(XW), nrow=nr, ncol=nr)
  ## old:
  ## XW <- as.matrix(expand.grid(X=obs[, "X"], W=obs[, "W"]))
  ## thetaXW <- matrix(theta(XW), nrow=nr, ncol=nr)
  ## new:
  XW <- as.matrix(rbind(expand.grid(X=Xq.index, W=1:nr),
                        cbind(X=1:nr, W=1:nr)))
  thetaXW <- theta(XW)
  
  if (!cleverCovTheta) {
    eps <- getEpsilon(this);
    g <- getG(this, tabulate=TRUE);
    mu <- getMu(this, tabulate=TRUE);
    sigma2 <- getSigma2(this);
    psi <- getPsi(this);

    obs <- getObs(this, tabulate=TRUE);
    gW <- g(obs[, "W"])
    ## old:
    ## gW <- t(matrix(gW, nrow=nr, ncol=nr))
    ## new:
    gW <- c(rep(gW, each=length(Xq.value)), gW)
    
    muW <- mu(obs[, "W"])
    ## old:
    ## muW <- t(matrix(muW, nrow=nr, ncol=nr))
    ## new:
    muW <- c(rep(muW, each=length(Xq.value)), muW)

    rm(g, mu, obs)

    obs <- getObs(this)
    ## old:
    ## eg <- expand.grid(X=obs[, "X"], W=obs[, "W"])
    ## OBSTAB <- cbind(X=fX(eg), fW(eg))
    ## X <- matrix(fX(obs), nrow=nr, ncol=nr)
    ## new:
    eg.value <- expand.grid(X=Xq.value, W=obs[, "W"])
    eg.index <- expand.grid(X=Xq.index, W=1:nr)
    OBSTAB.nondiag <- cbind(X=eg.value[, "X"], fW(eg.value))
    OBSTAB.diag <- cbind(X=fX(obs), W=fW(obs))
    OBSTAB <- rbind(OBSTAB.nondiag, OBSTAB.diag)
    X <- OBSTAB[, "X"]
    rm(obs)

    ## old:
    ## devThetaAuxXW <- matrix(devThetaAux(OBSTAB), nrow=nr, ncol=nr)
    ## devXW <- devThetaAuxXW * (X - muW/gW*(X==0))/sigma2;
    ## new:
    devThetaAuxXW <- devThetaAux(OBSTAB)
    devXW <- devThetaAuxXW * (X - muW/gW*(X==0))/sigma2
    
    if (!exact) { ## do not use clever covariate nor exact expression
      theta1XW <- thetaXW + eps*devXW
    } else { ## do not use clever covariate, but exact expression
      ## old:
      ## theta0W <- t(matrix(theta0(1:nr), nrow=nr, ncol=nr))
      ## new:
      theta0W <- theta0(1:nr)
      theta0W <- c(rep(theta0W, each=length(Xq.value)),
                   theta0W)
      term2 <- ( X*(thetaXW - theta0W -X*psi) )/sigma2;
      numerator <- thetaXW + eps*(devXW+term2*thetaXW);
      denominator <- 1 + eps*term2;
      theta1XW <- numerator/denominator;
    }
  } else { ## use clever covariate (hence, exact expression)
    if (!exact) {
      throw("Can't use a clever covariate for updating 'theta' if 'exact==FALSE' !");
    }
    eps <- getEpsilonTheta(this);
    dev <- getHTheta(this, tabulate=TRUE);
    ## old:
    ## devXW <- dev(expand.grid(X=1:nr, W=1:nr))
    ## devXW <- matrix(devXW, nrow=nr, ncol=nr)
    ## theta1XW <- thetaXW + eps*dev(XW)
    ## new:
    if (FALSE) {
      OBSTAB.nondiag <- expand.grid(X=Xq.index, W=1:nr)
      OBSTAB.diag <- cbind(X=1:nr, W=1:nr)
      OBSTAB <- rbind(OBSTAB.nondiag, OBSTAB.diag)
      theta1XW <- thetaXW + eps*dev(OBSTAB)
    } else {## should be the same
      theta1XW <- thetaXW + eps*dev(XW)
    }
  }
  ## old: 
  ## ## retrieving 'theta01W' from 'theta1XW'
  ## obs <- getObs(this, tabulate=FALSE)
  ## idx <- which(obs[, "X"]==0)[1]
  ## theta01W <- theta1XW[idx, ]
  ## new:
  whichXqIsZero <- which(eg.index[, "X"]==Xq.index[Xq.value==0])
  theta01W <- theta1XW[whichXqIsZero]

  ## old:
  ## THETA1TAB <- theta1XW
  ## new:
  ## making the sparse matrix
  THETA1TAB <- sparseMatrix(i=eg.index[, "X"], j=eg.index[, "W"],
                            x=as.vector(theta1XW[1:nrow(eg.index)]),
                            dims=c(nr, nr))
  THETA1TAB[cbind(1:nr, 1:nr)] <- as.vector(theta1XW[nrow(eg.index)+(1:nr)])
  theta1tab <- function(xiwj) {
    stopifnot(is.matrix(xiwj) && ncol(xiwj)==2 && is.integer(xiwj))
    THETA1TAB[xiwj]
  }
  THETA01TAB <- theta01W
  theta01tab <- function(wi) {
    ## stopifnot(is.matrix(wi) && ncol(wi)==1 && is.integer(wi))
    THETA01TAB[wi]
  }
  
  setThetaTab(this, theta1tab, theta01tab);
})

############################################################################
## HISTORY:
## 2014-02-07
## o Created.
############################################################################

