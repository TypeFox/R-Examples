#' Construct a posterior simulation object for the CTCRW state vectors
#' 

#' 
#' The \code{crwSimulator} function uses a fitted model object from
#' \code{crwMLE} and a set of prediction times to construct a list from which
#' \code{\link{crwPostIS}} will draw a sample from either the posterior
#' distribution of the state vectors conditional on fitted parameters or a full
#' posterior draw from an importance sample of the parameters.
#' 

#' 
#' The crwSimulator function produces a list and preprocesses the necessary
#' components for repeated track simulation from a fitted CTCRW model from
#' \code{\link{crwMLE}}. The \code{method} argument can be one of \code{"IS"}
#' or \code{"quadrature"}. If method="IS" is chosen standard importance
#' sampling will be used to calculate the appropriate weights via t proposal
#' with df degrees of freedom.  If df=Inf (default) then a multivariate normal
#' distribution is used to approximate the parameter posterior. If
#' \code{method="quadrature"}, then a regular grid over the posterior is used
#' to calculate the weights. The argument \code{grid.eps} controls the
#' quadrature grid. The arguments are approximately the upper and lower limit
#' in terms of standard deviations of the posterior.  The default is
#' \code{grid.eps}, in units of 1sd. If \code{object.crwFit} was fitted with
#' \code{crwArgoFilter}, then the returned list will also include \code{p.out},
#' which is the approximate probability that the observation is an outlier.
#' 
#' @param object.crwFit A model object from \code{\link{crwMLE}}.
#' @param predTime vector of additional prediction times.
#' @param method Method for obtaining weights for movement parameter samples
#' @param parIS Size of the parameter importance sample
#' @param df Degrees of freedom for the t approximation to the parameter
#' posterior
#' @param grid.eps Grid size for \code{method="quadrature"}
#' @param crit Criterion for deciding "significance" of quadrature points
#' (difference in log-likelihood)
#' @param scale Scale multiplier for the covariance matrix of the t
#' approximation
#' @return
#' 
#' List with the following elements:
#' 
#' \item{x}{Longitude coordinate with NA at prediction times}
#' 
#' \item{y}{Similar to above for latitude}
#' 
#' \item{locType}{Indicates prediction types with a "p" or observation times
#' with an "o"} \item{P1.y}{Initial state covariance for latitude}
#' 
#' \item{P1.x}{Initial state covariance for longitude}
#' 
#' \item{a1.y}{Initial latitude state}
#' 
#' \item{a1.x}{Initial longitude state}
#' 
#' \item{n.errX}{number of longitude error model parameters}
#' 
#' \item{n.errY}{number of latitude error model parameters}
#' 
#' \item{delta}{vector of time differences}
#' 
#' \item{driftMod}{Logical. indicates random drift model}
#' 
#' \item{stopMod}{Logical. Indicated stop model fitted}
#' 
#' \item{stop.mf}{stop model design matrix}
#' 
#' \item{err.mfX}{Longitude error model design matrix}
#' 
#' \item{err.mfY}{Latitude error model design matrix}
#' 
#' \item{mov.mf}{Movement model design matrix}
#' 
#' \item{fixPar}{Fixed values for parameters in model fitting}
#' 
#' \item{Cmat}{Covaraince matrix for parameter sampling distribution}
#' 
#' \item{Lmat}{Cholesky decomposition of Cmat}
#' 
#' \item{par}{fitted parameter values}
#' 
#' \item{N}{Total number of locations}
#' 
#' \item{loglik}{log likelihood of the fitted model}
#' 
#' \item{Time}{vector of observation times}
#' 
#' \item{coord}{names of coordinate vectors in original data}
#' 
#' \item{Time.name}{Name of the observation times vector in the original data}
#' 
#' \item{thetaSampList}{A list containing a data frame of parameter vectors and
#' their associated probabilities for a resample}
#' @author Devin S. Johnson
#' @seealso See \code{demo(northernFurSealDemo)} for example.
#' @export
crwSimulator = function(
  object.crwFit, 
  predTime=NULL, 
  method="IS", 
  parIS=1000, 
  df=Inf, 
  grid.eps=1, 
  crit=2.5, 
  scale=1) {
  ## Model definition/parameters ##
  data <- object.crwFit$data
  driftMod <- object.crwFit$random.drift
  mov.mf <- object.crwFit$mov.mf
  activity <- object.crwFit$activity
  err.mfX <- object.crwFit$err.mfX
  err.mfY <- object.crwFit$err.mfY
  rho = object.crwFit$rho
  par <- object.crwFit$par
  n.errX <- object.crwFit$n.errX
  n.errY <- object.crwFit$n.errY
  n.mov <- object.crwFit$n.mov
  tn <- object.crwFit$Time.name
  
  ## Data setup ##
  if (!is.null(predTime)) {
    if(min(predTime) <  data[1, tn]) {
      warning("Predictions times given before first observation!\nOnly those after first observation will be used.")
      predTime <- predTime[predTime>=data[1,tn]]
    }
    origTime <- data[, tn]
    if (is.null(data$locType)) data$locType <- "o"
    predData <- data.frame(predTime, "p")
    names(predData) <- c(tn, "locType")
    data <- merge(data, predData,
                  by=c(tn, "locType"), all=TRUE)
    dups <- duplicated(data[, tn]) #& data[,"locType"]==1
    data <- data[!dups, ]
    mov.mf <- as.matrix(expandPred(x=mov.mf, Time=origTime, predTime=predTime))
    if (!is.null(activity)) activity <- as.matrix(expandPred(x=activity, Time=origTime, predTime=predTime))
    if (!is.null(err.mfX)) err.mfX <- as.matrix(expandPred(x=err.mfX, Time=origTime, predTime=predTime))
    if (!is.null(err.mfY)) err.mfY <- as.matrix(expandPred(x=err.mfY, Time=origTime, predTime=predTime))
    if (!is.null(rho)) rho <- as.matrix(expandPred(x=rho, Time=origTime, predTime=predTime))
  }
  data$locType[data[,tn]%in%predTime] <- 'p'
  delta <- c(diff(data[, tn]), 1)
  a = object.crwFit$initial.state$a
  P = object.crwFit$initial.state$P
  y = as.matrix(data[,object.crwFit$coord])
  noObs <- as.numeric(is.na(y[,1]) | is.na(y[,2]))
  y[noObs==1,] = 0
  N = nrow(y)
  out <- list(y=y, noObs=noObs, P=P, a=a, n.errX=n.errX, n.errY=n.errY, n.mov=n.mov,
              delta=delta, driftMod=driftMod, activity=activity, err.mfX=err.mfX,
              err.mfY=err.mfY, mov.mf=mov.mf, rho=rho, fixPar=object.crwFit$fixPar,
              Cmat=object.crwFit$Cmat, locType=data$locType,
              par=object.crwFit$par, nms=object.crwFit$nms, N=nrow(data), lower=object.crwFit$lower, 
              upper=object.crwFit$upper,
              loglik=object.crwFit$loglik, Time=data[,tn], Time.name=tn, 
              coord=object.crwFit$coord, prior=object.crwFit$prior)
  class(out) <- 'crwSimulator'
  if(parIS>1 & object.crwFit$need.hess==TRUE) out <- crwSamplePar(out, method=method, size=parIS, df=df, grid.eps=grid.eps, crit=crit, scale=scale)
  return(out)
}
