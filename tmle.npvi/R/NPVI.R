setConstructorS3("NPVI", function(obs=matrix(nrow=0, ncol=3, dimnames=list(NULL, c("W", "X", "Y"))),
                                  f=identity, nMax=10L,
                                  gmin=0.01, gmax=1-gmin, mumin=-Inf,
                                  mumax=Inf, thetamin=-Inf, thetamax=Inf,
                                  family=c("parsimonious", "gaussian"), tabulate=TRUE,
                                  stoppingCriteria=list(mic=0.03, div=0.01, psi=0.1),
                                  conf.level=0.95,
                                  ...,
                                  verbose=FALSE) {
  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ## Validate arguments
  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ## Argument 'obs':
  obs <- validateArgumentObs(obs, allowIntegers=FALSE);
  
  ## Argument 'f':
  if (!((mode(f)=="function") && (f(0)==0))) {
    throw("Argument 'f' must be a function such that f(0)=0.")
  }

  ## Argument 'nMax':
  nMax <- Arguments$getInteger(nMax, c(10, Inf));
  
  ## Arguments 'gmin' and 'gmax':
  gmin <- Arguments$getNumeric(gmin);
  gmax <- Arguments$getNumeric(gmax);
  if (gmin>=1 | gmin<=0) {
    throw("Argument 'gmin' must belong to ]0,1[");
  }
  if (gmax>=1 | gmax<=0) {
    throw("Argument 'gmax' must belong to ]0,1[");
  }
  if (gmin>=gmax) {
    throw("Argument 'gmin' must be smaller than argument 'gmax'.")
  }

  ## Arguments 'mumin' and 'mumax':
  mumin <- Arguments$getNumeric(mumin);
  mumax <- Arguments$getNumeric(mumax);
  if (mumin>=mumax) {
    throw("Argument 'mumin' must be smaller than argument 'mumax'.")
  }

  ## Arguments 'thetamin' and 'thetamax':
  thetamin <- Arguments$getNumeric(thetamin);
  thetamax <- Arguments$getNumeric(thetamax);
  if (thetamin>=thetamax) {
    throw("Argument 'thetamin' must be smaller than argument 'thetamax'.")
  }

  ## Argument 'family'
  family <- match.arg(family)

  ## Argument 'tabulate'
  tabulate <- Arguments$getLogical(tabulate)

  ## Argument 'stoppingCriteria'
  mic.tol <- Arguments$getNumeric(stoppingCriteria$mic)
  attr(mic.tol, "label") <- "scaled empirical mean of estimating function"
  div.tol <- Arguments$getNumeric(stoppingCriteria$div)
  attr(div.tol, "label") <- "TV distance between P_n^k and P_n^{k+1}"
  psi.tol <- Arguments$getNumeric(stoppingCriteria$psi)
  attr(psi.tol, "label") <- "change between successive values of \"psi\""

  ## Argument 'conf.level'
  conf.level <- Arguments$getNumeric(conf.level, c(0, 1))

  ## Arguments '...':
  args <- list(...);
  if (length(args) > 0) {
    argsStr <- paste(names(args), collapse=", ");
    throw("Unknown arguments: ", argsStr);
  } 

  if (tabulate & family=="gaussian") {
    throw("Unauthorized, because cannot tabulate all functions if family is 'gaussian'.")
  }
  if (!tabulate & family=="parsimonious") {
    throw("Unauthorized, because not implemented (it is a sub-optimal combination of options).")
  }

  obs[, "X"] <- X <- f(obs[, "X"])
  if (length(X)==0) {
    Xq <- data.frame(value=numeric(0), index=integer(0))
  } else {
    Xneq0 <- X[X!=0]
    Xq <- unique(quantile(Xneq0, type=1, probs=seq(0, 1, length=nMax-1)))
    if (length(setdiff(Xq, Xneq0))) {
      throw("In 'NPVI': components of 'Xq' must be observed values of 'X'...")
    }
    Xq.idx <- match(Xq, X)
    Xq0.idx <- which(X==0)[1]
    Xq <- data.frame(value=c(0, Xq), index=c(Xq0.idx, Xq.idx))
  }

  Y <- obs[, "Y"]
  if (length(Y)==0) {
    Yq <- data.frame(value=numeric(0), index=integer(0))
  } else {
    Yq <- unique(quantile(Y, type=1, probs=seq(0, 1, length=nMax)))
    if (length(setdiff(Yq, Y))) {
      throw("In 'NPVI': components of 'Yq' must be observed values of 'Y'...")
    }
    Yq.idx <- match(Yq, Y)
    Yq <- data.frame(value=Yq, index=Yq.idx)
  }

  
  theW <- setdiff(colnames(obs), c("X", "Y"))
  if (!tabulate & length(theW)>1) {
    throw("Multivariate 'W' handled only if 'tabulate' is TRUE")
  }
  
  nms <-  c("eps", "lli", "mic1", "epsT", "lliT", "mic2", "psi", "psi.sd", "psiPn", "psiPn.sd", "mic", "div", "sic", "phi", "sicAlt")
  history <- matrix(NA, 0, length(nms));
  colnames(history) <- nms

  conv <- NA
  attr(conv, "msg") <- character(0)
  
  extend(Object(), "NPVI",
         .obs=obs, #.flavor=flavor,
         .Xq=Xq, .Yq=Yq,
         .g=NULL, .mu=NULL, .muAux=NULL, .theta=NULL, .theta0=NULL, .weightsW=rep(1, nrow(obs)), 
         .gtab=NULL, .mutab=NULL, .muAuxtab=NULL, .thetatab=NULL, .theta0tab=NULL,
         .sigma2=NA, .psi=NA, .psi.sd=NA, .psiPn=NA, .psiPn.sd=NA,
         .gmin=gmin, .gmax=gmax, .mumin=mumin, .mumax=mumax,
         .thetamin=thetamin, .thetamax=thetamax,
         .family=family, .tabulate=tabulate, .epsilon=NA,
         .epsilonTheta=NA, .logLikIncr=NA, .logLikIncrTheta=NA, .div=NA,
         .efficientInfluenceCurve=matrix(NA, 0, 3), .history=history, .step=0,
         .stoppingCriteria=list(mic=mic.tol, div=div.tol, psi=psi.tol),
         .conv=conv, .conf.level=conf.level
         );
})

setMethodS3("getXq", "NPVI", function(this, ...) {
  this$.Xq;
})

setMethodS3("getYq", "NPVI", function(this, ...) {
  this$.Yq;
})


## setMethodS3("getFlavor", "NPVI", function(this, ...) {
##   this$.flavor;
## })

setMethodS3("getMicTol", "NPVI", function(this, ...) {
  this$.stoppingCriteria$mic;
})

setMethodS3("getDivTol", "NPVI", function(this, ...) {
  this$.stoppingCriteria$div;
})

setMethodS3("getPsiTol", "NPVI", function(this, ...) {
  this$.stoppingCriteria$psi;
})

setMethodS3("getConfLevel", "NPVI", function(this, ...) {
  this$.conf.level;
})

setMethodS3("getHistory", "NPVI", function(#Returns History of TMLE Procedure
  ### Returns the 'history' of the TMLE procedure.
  this,
### An object of class \code{TMLE.NPVI}.
  ...
### Not used.
  ) {
  ##alias<< getHistory
  ##seealso<< tmle.npvi
  this$.history;
  ###  Returns a \code{numeric}  \code{matrix} which encapsulates a summary of
  ###     the  TMLE procedure.  If \eqn{k} successive  updates were performed,
  ###       then   the   \code{matrix}   has   either   \eqn{k+1}   rows   (if
  ###      \code{cleverCovTheta}  was  set  to  \code{FALSE} in  the  call  to
  ###       \code{tmle.npvi})   or    \code{2k+1}   rows   (otherwise).    The
  ###     \code{matrix} has 14 columns:
  ### \itemize{
  ###       \item{\code{"eps"},   values of  the unique
  ###        fluctuation  parameter   (if  \code{cleverCovTheta}  was  set  to
  ###        \code{FALSE} in  the call  to \code{tmle.npvi}),  or  
  ###        values of the  parameter involved in the fluctuation of
  ###        the   joint  distribution  of  \eqn{(X,W)}   during  each  update
  ###       (otherwise).  }
  ###     \item{\code{"lli"},     increases in likelihood
  ###      yielded  by  each  update  (if  \code{cleverCovTheta}  was  set  to
  ###      \code{FALSE}  in the  call  to  \code{tmle.npvi}),  or 
  ###      increases in likelihood yielded by the fluctuation of the
  ###     joint distribution of \eqn{(X,W)} during each update (otherwise).}
  ###        \item{\code{"mic1"},   empirical means  of the  first
  ###     component of the efficient influence  curve at each step of the TMLE
  ###     procedure.}
  ### \item{\code{"epsT"},   values of  the fluctuation
  ### parameter involved in the fluctuation of the conditional distribution of
  ### \eqn{Y}  given \eqn{(X,W)} during each  update (if \code{cleverCovTheta}
  ### was  set to \code{TRUE} in  the call to  \code{tmle.npvi}), or \code{NA}
  ### (otherwise).}
  ###     \item{\code{"lliT"},  successive increases in likelihood
  ###       yielded by  the  fluctuation of  the  conditional distribution  of
  ###         \eqn{Y}    given   \eqn{(X,W)}    during    each   update    (if
  ###       \code{cleverCovTheta}  was  set  to  \code{TRUE} in  the  call  to
  ###      \code{tmle.npvi}), or \code{NA} (otherwise).}
  ### \item{\code{"mic2"},  empirical means  of the  second
  ###     component of the efficient influence  curve at each step of the TMLE
  ###     procedure.}
  ###   \item{\code{"psi"},   increasingly targeted
  ###     estimators \eqn{\Psi(P_n^k)} of  the  parameter  of  interest. The last one is the TMLE.  Their  computation  
  ###       involves  simulation  of  \code{B}  iid  copies  of
  ###    \eqn{(X,W)} under \eqn{P_n^k}. }
  ### \item{\code{"psi.sd"},  estimated standard deviations of
  ###  the   increasingly targeted  estimators of  the  parameter of
  ###  interest. The last one corresponds to the TMLE. The  computation involves  the  same \code{B}  iid copies  of
  ### \eqn{(X,W)} as above.}
  ###  \item{\code{"psiPn"},  same as \code{"psi"} except  that the *observed*
  ###   \eqn{(X_i,W_i)}  are  used  instead  of simulated  copies  drawn  from
  ###  \eqn{P_n^k}. Of course, \code{"psi"} must be favored.}
  ###   \item{\code{"psiPn.sd"},  same  as  \code{"psi.sd"}  except  that  the
  ###  *observed*  \eqn{(X_i,W_i)} are used instead of  simulated copies drawn
  ###  from \eqn{P_n^k}. Of course, \code{"psi.sd"} must be favored.}
  ###  \item{\code{"mic"},   empirical  means  of  the  efficient
  ### influence curve at each step of the TMLE procedure. This column is the sum of the \code{"mic1"} and \code{"mic2"} columns.}
  ### \item{\code{"div"},  total variation  distances between
  ### each pair  of successive distributions constructed in  the course of the
  ### TMLE procedure. }
  ### \item{\code{"sic"},  estimated standard deviations   of  the  efficient
  ### influence curve at each step of the TMLE procedure.}
  ###\item{\code{"phi"},    non-parametric     substitution    estimator    of
  ###\eqn{\phi=\Phi(P)}             where            \deqn{\Phi(P)            =
  ###\frac{E_P[f(X)Y]}{E_P[f(X)^2]},}{\Phi(P)  =  E_P[f(X)Y]  /  E_P[f(X)^2],}
  ###with  \eqn{P}  the  distribution   of  the  random  vector  \eqn{(W,X,Y)}. The alternative parameter \eqn{\phi} should be interpreted as the counterpart of \eqn{\psi} which neglects \eqn{W}. }
  ### \item{\code{"sicAlt"},  estimated standard deviations   of  the  efficient
  ### influence curve of \eqn{\Psi - \Phi} at each step of the TMLE procedure.}
###   }
})

setMethodS3("updateHistory", "NPVI", function(this, ...) {
  history <- tmle.npvi::getHistory.NPVI(this);
  psi <- getPsi(this);
  psi.sd <- getPsiSd(this);  
  epsilon <- getEpsilon(this);
  logLikIncr <- getLogLikIncr(this);
  epsilonTheta <- getEpsilonTheta(this);
  logLikIncrTheta <- getLogLikIncrTheta(this);
  eic <- getEfficientInfluenceCurve(this);
  mic <- apply(eic, 2, mean);
  sic <- getSic(this);
  div <- getDivergence(this)
  psiPn <- getPsiPn(this);
  psiPn.sd <- getPsiPnSd(this);

  phi <- getPhi(this)
  sicAlt <- getSicAlt(this)
  
  currStep <- c(epsilon, logLikIncr, mic[1], epsilonTheta, logLikIncrTheta, mic[2], psi, psi.sd, psiPn, psiPn.sd, mic[3], div, sic, phi, sicAlt);
  step <- getStep(this);
  rownames.history <- rownames(history);
  history <- rbind(history, matrix(currStep, 1, length(currStep)));
  rownames(history) <- c(rownames.history,
                         paste("step", as.character(step), sep=""));
  
  this$.history <- history;
})

setMethodS3("getPsi", "NPVI", function(#Returns Current Estimator
### Returns the current value of the estimator.
    this,
### An object of class \code{TMLE.NPVI}.
    ...
### Not used.
    ) {
  ##alias<< getPsi
  ##seealso<< tmle.npvi, getHistory, getPsiSd
  this$.psi;
  ### Retrieves  the current value  of the estimator \eqn{\Psi(P_n^k)}  of the
  ### parameter  of interest. Its computation involves  simulation of a large number of 
  ### iid copies of \eqn{(X,W)} under \eqn{P_n^k}.
})

setMethodS3("getPsiSd", "NPVI", function(#Returns Current Estimated Standard Deviation of the Estimator 
### Returns the current value of the estimated standard deviation of the current estimator. 
    this,
### An object of class \code{TMLE.NPVI}.
    ...
### Not used.
    ) {
  ##alias<< getPsiSd
  ##seealso<< tmle.npvi, getHistory, getPsi
  this$.psi.sd;
  ### Retrieves  the estimated standard deviation of the current estimator \eqn{\Psi(P_n^k)}  of the
  ### parameter  of interest. Its computation involves  simulation of a large number of 
  ### iid copies of \eqn{(X,W)} under \eqn{P_n^k}.  
})

setMethodS3("getPhi", "NPVI", function(this, ...) {
  obs <- getObs(this)
  fX <- getFX(this)
  fY <- getFY(this)
  X <- fX(obs)
  Y <- fY(obs)
  sX2 <- mean(X^2)
  ## phi: estimator of phi_0 
  mean(X*Y)/sX2
})

setMethodS3("getSic", "NPVI", function(this, ...) {
  eic <- getEfficientInfluenceCurve(this);
  mic <- apply(eic, 2, mean);
  sd(eic[, 3])
})

setMethodS3("getSicAlt", "NPVI", function(this, ...) {
  obs <- getObs(this)
  fX <- getFX(this)
  fY <- getFY(this)
  X <- fX(obs)
  Y <- fY(obs)
  eic <- getEfficientInfluenceCurve(this)
  ## phi: estimator of phi_0 
  phi <- getPhi(this)
  ## sicAlt: estimated standard deviation to perform test of "psi_0 = phi_0"
  sX2 <- mean(X^2)
  infCurvePhi <- (X*Y-phi*X^2)/sX2
  sicAlt <- sd(eic[, 3]-infCurvePhi)
  sicAlt
})


setMethodS3("getPsiPn", "NPVI", function(this, name, ...) {
  this$.psiPn;
})

setMethodS3("getPsiPnSd", "NPVI", function(this, name, ...) {
  this$.psiPn.sd;
})

setMethodS3("getStep", "NPVI", function(this, name, ...) {
  this$.step;
})

setMethodS3("getDivergence", "NPVI", function(this, name, ...) {
  this$.div;
})

setMethodS3("setDivergence", "NPVI", function(this, div, ...) {
  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ## Validate arguments
  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ## Argument 'div':
  div <- Arguments$getNumeric(div)
  
  this$.div <- div;
})

setMethodS3("setConfLevel", "NPVI", function(#Sets Confidence Level
### Sets the confidence level of a \code{TMLE.NPVI} object.
    this,
### An object of class \code{TMLE.NPVI}.
    confLevel,
### A \code{numeric}, confidence interval level.
    ...
### Not used.
    ) {
  ##alias<< setConfLevel
  ##seealso<< as.character.NPVI
  
  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ## Validate arguments
  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ## Argument 'confLevel':
  conf.level <- Arguments$getNumeric(confLevel, range=c(0, 1))
  
  this$.conf.level <- conf.level;
})

setMethodS3("getFW", "NPVI", function(this, tabulate, ...) {
  if (missing(tabulate)) {
    tabulate <- getTabulate(this)
  }
  if (!tabulate) {
    ## 'W' is necessarily unidimensional
    fW <- function(obs) {
      obs[, "W", drop=FALSE];
    }
  } else {
    ## 'W' needs not be unidimensional
    obsT <- getObs(this, tabulate=FALSE);
    fW <- function(obs) {
      ## obsT[obs[, "W"], "W", drop=FALSE];
      extractW(obsT[obs[, "W"], ])
    }
  }
  fW
})


setMethodS3("getFX", "NPVI", function(this, tabulate, ...) {
  if (missing(tabulate)) {
    tabulate <- getTabulate(this)
  }
  if (!tabulate) {
    fX <- function(obs) {
      obs[, "X"];
    }
  } else {
    obsT <- getObs(this, tabulate=FALSE);
    fX <- function(obs) {
      obsT[obs[, "X"], "X"];
    }
  }
  fX
})

setMethodS3("getFY", "NPVI", function(this, tabulate, ...) {
  if (missing(tabulate)) {
    tabulate <- getTabulate(this)
  }
  if (!tabulate) {
    fY <- function(obs) {
      obs[, "Y"];
    }
  } else {
    obsT <- getObs(this, tabulate=FALSE);
    fY <- function(obs) {
      obsT[obs[, "Y"], "Y"];
    }
  }
  fY
})

setMethodS3("getObs", "NPVI", function(#Retrieves the Observations
### Retrieves the \code{matrix} of observations involved in the TMLE procedure.
    this,
### An object of class \code{TMLE.NPVI}.
    tabulate,
### A \code{logical}, to  specify whether it is the original  data set that is
### retrieved (if \code{FALSE}) or a  tabulated version of it (otherwise), for
### internal use only.  If \code{tabulate} is missing then  the value attached
### to the input object is used.
    ...
### Not used.
    ) {
  ##alias<< getObs
  ##seealso<< tmle.npvi
  if (missing(tabulate)) {
    tabulate <- getTabulate(this);
  }
  if (!tabulate) {
    obs <- this$.obs;
  } else {
    n <- nrow(this$.obs)
    obs <- matrix(1:n, nrow=n, ncol=3);
    colnames(obs) <- c("W", "X", "Y");
  }
  obs
### Either the original data set involved in the TMLE procedure or a tabulated
### version of it.
})

setMethodS3("getFamily", "NPVI", function(this, ...) {
  this$.family;
})

setMethodS3("setFamily", "NPVI", function(this, family=c("parsimonious", "gaussian"), ...) {
  ## Argument
  family <- match.arg(family);
  this$.family <- family;
})

setMethodS3("getTabulate", "NPVI", function(this, ...) {
  this$.tabulate;
})

setMethodS3("getHTheta", "NPVI", function(this, tabulate, ...) {
  if (missing(tabulate)) {
    tabulate <- getTabulate(this)
  }
  mu <- getMu(this, tabulate)
  g <- getG(this, tabulate)
  sigma2 <- getSigma2(this)
  fX <- getFX(this, tabulate);
    
  HTheta <- function(XW) {
    X <- fX(cbind(XW, Y=NA))
    W <- XW[, "W"]
    (X - (X==0)*mu(W)/g(W))/sigma2
  }
  HTheta;
})

setMethodS3("getSigma2", "NPVI", function(this, ...) {
  this$.sigma2;
})

setMethodS3("setSigma2", "NPVI", function(this, sigma2, ...) {
  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ## Validate arguments
  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ## Argument 'sigma2':
  sigma2 <- Arguments$getNumeric(sigma2)

  this$.sigma2 <- sigma2;
})

setMethodS3("as.character", "NPVI", function(#Returns a Description
### Returns a short string describing the NPVI object.
    x,
### An object of class \code{TMLE.NPVI}.
    ...
### Not used.
    ) {
  ##alias<< as.character
  ##seealso<< tmle.npvi
  this <- x;  ## To please R CMD check
  s <- sprintf("%s object:", class(this)[1]);
  s <- c(s, "")

  ##   switched to from 'superLearning' to 'learning'?
  flag <- attr(this, "flag")
  if (!is.null(flag)) {
    s <- c(s, flag, "")
  } 

  
  ## sample size
  s <- c(s, sprintf("Sample size: %s", nrow(getObs(this))))
  s <- c(s, "")
  
  ## psi
  s <- c(s, sprintf("Estimator of psi:\t\t%s", signif(getPsi(this), 3)));

  ## std error
  sic <- getSic(this)
  s <- c(s, sprintf("Estimated standard error:\t%s", signif(sic, 3)))
  s <- c(s, "")
  
  ## number of iterations
  step <- getStep(this)
  
  ## convergence ?
  conv <- getConv(this)
  if (!is.na(conv)) {
    if (conv) {
      msg <- paste("Convergence reached after", step, "iteration(s) because:")
      msg <- paste("\n", msg, "\n\t",
                   paste(attr(conv, "msg"), collapse="\n\t"),
                   sep="")
    } else {
      msg <- paste("Convergence not reached after", step, "iteration(s)")
    }

    mic <- getMicTol(this)
    micLab <- attr(mic, "label")
    div <- getDivTol(this)
    divLab <- attr(div, "label")
    psi <- getPsiTol(this)
    psiLab <- attr(psi, "label")
    
    msg2 <- paste("Convergence criteria: ",
                  "\n- ", micLab, "\t\t< ", mic, 
                  "\n- ", divLab, "\t\t< ", div, 
                  "\n- ", psiLab, "\t\t< ", psi, sep="")
    rm(psi, mic, div)
    s <- c(s, msg2, msg)
    s <- c(s, "")
  }
  
  ## confidence intervals
  psi <- getPsi(this)
  alpha <- 1-getConfLevel(this)
  n <- nrow(getObs(this))
  CI <- psi+c(-1, 1)*sic*qnorm(1-alpha/2)/sqrt(n)
  CI <- signif(CI, 3)
  s <- c(s, sprintf("%s-confidence interval:\t[%s, %s]", 1-alpha, CI[1], CI[2]))
  
  ## tests
  ts1 <- sqrt(n)*(psi-0)/sic
  pval1 <- 2*(1-pnorm(abs(ts1)))
  s <- c(s, sprintf("Test of \"psi(P_0)=0\":\t\tp-value = %s", signif(pval1, 3)))

  phi <- getPhi(this)
  ts2 <- sqrt(n)*(psi-phi)/getSicAlt(this)
  pval2 <- 2*(1-pnorm(abs(ts2)))
  s <- c(s, sprintf("Test of \"psi(P_0)=phi(P_0)\":\tp-value = %s",
                    signif(pval2, 3)),
         sprintf(" (estimated phi(P_0):\t%s)", signif(phi, 3)))

  class(s) <- "GenericSummary";
  s;
### A character string summarizing the content of the object. The summary contains:
### \itemize{
### \item{The sample size of the data set involved in the TMLE procedure.}
### \item{The value of the TMLE and its estimated standard error.}
### \item{A reminder of  the tuning of the stopping criteria,  and a report on
### the convergence of the TMLE procedure (see \code{\link{tmle.npvi}}).}
### \item{A confidence interval with default level of 95% (the level can be changed by using \code{\link{setConfLevel}}).}
### \item{The \eqn{p}-value of the two-sided test of ``\eqn{\Psi(P_0)=0}''.}
### \item{The \eqn{p}-value of the two-sided test of ``\eqn{\Psi(P_0)=\Phi(P_0)}'', with the estimated value of \eqn{\Phi(P_0)}.}
### }
}, private=TRUE)

setMethodS3("updateSigma2", "NPVI", function(this, dev, ...) {
  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ## Validate arguments
  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ## Argument 'dev':
  dev <- Arguments$getNumeric(dev);

  eps <- getEpsilon(this);  
  sigma2 <- getSigma2(this);
  
  sigma21 <- sigma2 + eps*dev;
  ## sigma21 is positive because 'eps' is upper bounded by 1/supremum of
  ## absolute value of efficient influence curve
  setSigma2(this, sigma21)
})

setMethodS3("updateEfficientInfluenceCurve", "NPVI", function(this, ..., verbose=FALSE) {
  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ## Retrieve arguments
  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  obs <- getObs(this)
  g <- getG(this);
  mu <- getMu(this);
  theta <- getTheta(this);
  theta0 <- getTheta0(this);
  fX <- getFX(this)
  fY <- getFY(this)
  sigma2 <- getSigma2(this);
  psi <- getPsi(this);
  
  thetaXW <- theta(obs[, c("X", "W")]);
  W <- obs[, "W", drop=FALSE]
  theta0W <- theta0(W);
  muW <- mu(W);
  gW <- g(W);

  X <- fX(obs)
  Y <- fY(obs)
          
  D1 <- X * (thetaXW - theta0W - X * psi);
  D2 <- (Y - thetaXW) * (X - muW/gW*(X==0));
  verbose && summary(verbose, D1);
  verbose && summary(verbose, D2);
  verbose && print(verbose, sigma2);

  eic1 <- D1 / sigma2;
  eic2 <- D2 / sigma2;
  eic <- eic1 + eic2;
  verbose && summary(verbose, eic);

  this$.efficientInfluenceCurve <- cbind(eic1, eic2, eic);
})

setMethodS3("estimateEpsilon", "NPVI", function(this, cleverCovTheta, bound=1e-1, ..., verbose=FALSE) {
  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ## Validate arguments
  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  ## Argument 'cleverCovTheta'
  cleverCovTheta <- Arguments$getLogical(cleverCovTheta);

  eic <- getEfficientInfluenceCurve(this, verbose=verbose);
  if (cleverCovTheta) {
    eic <- eic[, "eic1"]
  } else {
    eic <- eic[, "eic"];
  }
  verbose && summary(verbose, eic);
  
  if (sum(abs(eic)==Inf, na.rm=TRUE)>0) {
    throw("Infinite values in estimated efficient influence curve");
    ## eic[abs(eic)==Inf] <- NA;
  }

  theMin <- min(eic)
  theMax <- max(eic)
  if (theMin > 0) {
    interval <- c(-1/(1.001*theMax), 1e3)
  }
  if (theMax < 0) {
    interval <- c(-1e3, -1/(1.001*theMin))
  }
  if (theMin<=0 & theMax>=0) {
    interval <- c(-1/(1.001*theMax), -1/(1.001*theMin))
  }
  
  logLik <- function(epsilon) {
    sum(log(1 + epsilon * eic));
  }

  interval <- pmin(bound, pmax(-bound, interval))
  verbose && cat(verbose, "Optimization interval");
  verbose && print(verbose, interval);
  opt <- optimize(logLik, interval=interval, maximum=TRUE);

  names(opt) <- c("eps", "feps");
  eps <- opt$eps;
  attr(eps, "feps") <- opt$feps;

  eps;
})

setMethodS3("estimateEpsilonTheta", "NPVI", function(this, ..., verbose=FALSE) {
  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ## Validate arguments
  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ## Argument 'obs':
  obs <- getObs(this);
  fY <- getFY(this)
  
  theta <- getTheta(this);
  H <- getHTheta(this);

  XW <- obs[, c("X", "W")]
  HXW <- H(XW)
  residuals <- fY(obs)-theta(XW)
  
  eps <- mean(residuals*HXW)/mean(HXW^2);
  feps <- eps*sum(residuals*HXW)  - (eps^2) * sum(HXW^2)/2
  attr(eps, "feps") <- feps
  
  eps
})

setMethodS3("updateEpsilon", "NPVI", function(this, ..., verbose=FALSE) {
  eps <- estimateEpsilon(this, ..., verbose=verbose)

  this$.logLikIncr <- attr(eps, 'feps')
  attr(eps, 'feps') <- NULL
  this$.epsilon <- eps;
})

setMethodS3("updateEpsilonTheta", "NPVI", function(this, ..., verbose=FALSE) {
  epsTheta <- estimateEpsilonTheta(this, ..., verbose=verbose)

  this$.logLikIncrTheta <- attr(epsTheta, 'feps')
  attr(epsTheta, 'feps') <- NULL
  this$.epsilonTheta <- epsTheta;
})

setMethodS3("getEpsilon", "NPVI", function(this, ...) {
  this$.epsilon
})

setMethodS3("getEpsilonTheta", "NPVI", function(this, ...) {
  this$.epsilonTheta
})

setMethodS3("getLogLikIncr", "NPVI", function(this, ...) {
  this$.logLikIncr
})

setMethodS3("getLogLikIncrTheta", "NPVI", function(this, ...) {
  this$.logLikIncrTheta
})

setMethodS3("getEfficientInfluenceCurve", "NPVI", function(this, ...) {
  this$.efficientInfluenceCurve;
})

setMethodS3("getWeightsW", "NPVI", function(this, ...) {
  this$.weightsW
})

setMethodS3("setWeightsW", "NPVI", function(this, weightsW, ...) {
  ## Argument 'weightsW':
  weightsW <- Arguments$getNumerics(weightsW, range=c(0, Inf))
  nr <- nrow(getObs(this))
  if (length(weightsW) != nr) {
    throw("Length of 'weightsW' must match number of observations!")
  }
  this$.weightsW <- weightsW
})

setMethodS3("updateWeightsW", "NPVI", function(this, effICW, ...) {
  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ## Validate arguments
  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ## Argument 'effICW':
  if (!is.null(effICW) & (mode(effICW)!="function")) {
    throw("Argument \var{effICW} should be of mode 'function', not ", mode(effICW));
  }

  if (!is.null(effICW)) {
    obs <- getObs(this)
    eps <- getEpsilon(this);  
    weightsW <- getWeightsW(this)*(1+eps*effICW(obs[, "W"]))

    setWeightsW(this, weightsW)
  }
})

setMethodS3("getConv", "NPVI", function(this, ...) {
  this$.conv;
})

setMethodS3("updateConv", "NPVI", function(x, B, ...) {
  this <- x;  ## To please R CMD check

  n <- nrow(getObs(this))
  mic.tol <- getMicTol(this)
  div.tol <- getDivTol(this)
  psi.tol <- getPsiTol(this)
  ## extracting the relevant components of history
  ## (whether 'cleverCovTheta' is TRUE or FALSE)
  hist <- tmle.npvi::getHistory.NPVI(this)
  step <- as.integer(gsub("step([0-9]+)", "\\1", rownames(hist)))
  idxs <- c(which(diff(step)==1), length(step))
  hist <- tail(hist[idxs, ], 2)

  conv <- FALSE
  msg <- NULL
  
  ## Assessing values of scaled 'mic'
  if (!is.na(mic.tol)) {
    mic <- hist[2, "mic"]
    Zmic <- hist[2, "psi.sd"]/sqrt(n)
    scaledMic <- mic/Zmic
    micLab <- attr(mic.tol, "label")
    if (abs(scaledMic) < mic.tol) {
      conv <- TRUE;
      msg1 <- paste(micLab, " is within ", mic.tol, "-tolerance", sep="");
      msg <- c(msg, msg1)
    }
  }
  
  ## Testing value of 'div'
  if (!is.na(div.tol)) {
    div <- hist[2, "div"]
    divLab <- attr(div.tol, "label")
    if (!is.na(div)) {
      if (div < div.tol) {
        conv <- TRUE;
        msg1 <- paste(divLab, " is within ", div.tol, "-tolerance", sep="");
        msg <- c(msg, msg1)
      }
    }
  }

  ## Comparing successive values of 'psi'
  if (!is.na(psi.tol)) {
    psiLab <- attr(psi.tol, "label")
    psi <- hist[, "psi"]
    psi.sd <- hist[, "psi.sd"]
    crit <- abs(diff(psi))/sqrt(sum(psi.sd^2))*sqrt(B)
    if (crit < psi.tol) {
      conv <- TRUE;
      msg1 <- paste(psiLab, " are within ", psi.tol, "-tolerance", sep="");
      msg <- c(msg, msg1)
    }
  }

  if (!conv) {
    msg <- "Convergence not reached (yet)"
  }

  
  attr(conv, "msg") <- msg
  this$.conv <- conv
})

setMethodS3("getPValue", "NPVI", function(# Calculates a p-value from a NPVI object
### Calculates a p-value from a NPVI object
    this,
### An object of class \code{TMLE.NPVI}.
    wrt.phi=TRUE,
### A  \code{logical}  equal  to  \code{TRUE}  by default,  which  means  that
### \eqn{psi_n}  is  compared  with  \eqn{phi_n}.  Otherwise,  \eqn{psi_n}  is
### compared with 0.
    ...
### Not used.

){
    ##seealso<< tmle.npvi, getHistory, as.character.NPVI, getPValue.matrix
    nobs <- nrow(getObs(this))
    history <- getHistory.NPVI(this)
    getPValue(history, wrt.phi=wrt.phi, nobs=nobs)
### Returns the p-value of the two-sided test of
### ``\eqn{Psi(P_0)=Phi(P_0)}'' or ``\eqn{Psi(P_0)=0}'', according to
### the value of \code{wrt.phi}.
})

setMethodS3("getPValue", "matrix", function(# Calculates a p-value from a matrix object of type 'history'
### Calculates a p-value from a matrix object of type 'history'
    this,
### The \code{history} of a TMLE procedure.
    wrt.phi=TRUE,
### A  \code{logical}  equal  to  \code{TRUE}  by default,  which  means  that
### \eqn{psi_n}  is  compared  with  \eqn{phi_n}.  Otherwise,  \eqn{psi_n}  is
### compared with 0.
    nobs,
### An \code{integer}, the associated number of observations.
    ...
### Not used.
){
    ##alias<< getPValue
    ##seealso<< tmle.npvi, getHistory.NPVI, as.character.NPVI
    y <- this[nrow(this), ]
    psi <- y["psi"]
    if (wrt.phi) {
        phi <- y["phi"]
        se <- y["sicAlt"]/sqrt(nobs)
    } else {
        phi <- 0
        se <- y["psi.sd"]/sqrt(nobs)
    }
    pval <- 2*pnorm(abs(psi-phi), sd=se, lower.tail=FALSE)
    names(pval) <- "p.value"
    return(pval)
### Returns the p-value of the two-sided test of
### ``\eqn{Psi(P_0)=Phi(P_0)}'' of ``\eqn{Psi(P_0)=0}'', according to
### the value of \code{wrt.phi}.
})


############################################################################
## HISTORY:
## 2014-02-07
## o Created.
############################################################################

