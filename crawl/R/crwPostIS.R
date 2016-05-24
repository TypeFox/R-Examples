#' Simulate a value from the posterior distribution of a CTCRW model
#' 

#' 
#' The crwPostIS draws a set of states from the posterior distribution of a
#' fitted CTCRW model. The draw is either conditioned on the fitted parameter
#' values or "full" posterior draw with approximated parameter posterior
#' 

#' 
#' The crwPostIS draws a posterior sample of the track state matrices. If
#' fullPost was set to TRUE when the object.sim was build in
#' \link{crwSimulator} then a psuedo-posterior draw will be made by first
#' sampling a parameter value from a multivariate t distribution which
#' approximates the marginal posterior distribution of the parameters. The
#' covariance matrix from the fitted model object is used to scale the MVt
#' approximation. In addition, the factor "scale" can be used to further adjust
#' the approximation. Further, the parameter simulations are centered on the
#' fitted values.
#' 
#' To correct for the MVt approximation, the importance sampling weight is also
#' supplied. When calulating averages of track functions for Bayes estimates
#' one should use the importance sampling weights to calculate a weighted
#' average (normalizing first, so the weights sum to 1).
#' 
#' @param object.sim A crwSimulator object from \code{\link{crwSimulator}}.
#' @param fullPost logical. Draw parameter values as well to simulate full
#' posterior
#' @param df degrees of freedom for multivariate t distribution approximation
#' to parameter posterior
#' @param scale Extra scaling factor for t distribution approximation
#' @param thetaSamp If multiple parameter samples are available in object.sim,
#' setting \code{thetaSamp=n} will use the nth sample. Defaults to the last.
#' @return
#' 
#' List with the following elements:
#' 
#' \item{alpha.sim.y}{A matrix a simulated latitude state values}
#' 
#' \item{alpha.sim.x}{Matrix of simulated longitude state values}
#' 
#' \item{locType}{Indicates prediction types with a "p" or observation times
#' with an "o"} \item{Time}{Initial state covariance for latitude}
#' 
#' \item{loglik}{log likelihood of simulated parameter}
#' 
#' \item{par}{Simulated parameter value}
#' 
#' \item{log.isw}{non normalized log importance sampling weight}
#' @author Devin S. Johnson
#' @seealso See \code{demo(northernFurSealDemo)} for example.
#' @export

`crwPostIS` <-
  function(object.sim, fullPost=TRUE, df=Inf, scale=1, thetaSamp=NULL)
    ################################################################################
################################################################################
{
  if(!inherits(object.sim, 'crwSimulator')) stop("Argument needs to be of class 'crwSimulator'\nUse 'crwSimulator( )' to create")
  fixPar <- object.sim$fixPar
  Cmat <- object.sim$Cmat[is.na(fixPar),is.na(fixPar)]
  se <- sqrt(diag(Cmat))
  err.mfX <- object.sim$err.mfX
  err.mfY <- object.sim$err.mfY
  parMLE <- object.sim$par
  n2ll.mode <- -2*object.sim$loglik
  activity <- object.sim$activity
  driftMod <- object.sim$driftMod
  mov.mf <- object.sim$mov.mf
  y <- object.sim$y
  noObs = object.sim$noObs
  delta <- object.sim$delta
  a <- object.sim$a
  P <- object.sim$P
  n.errX <- object.sim$n.errX
  n.errY <- object.sim$n.errY
  rho = object.sim$rho
  n.mov <- object.sim$n.mov
  N <- object.sim$N
  lower <- object.sim$lower
  upper <- object.sim$upper 
  prior <- object.sim$prior
  eInd <- is.na(fixPar)
  ###
  ### Sample parameter vector
  ###
  if(fullPost & is.null(object.sim$thetaSampList)) {
    eps <- rmvtt(mu=rep(0,sum(eInd)), Sigma=scale*Cmat, df=df, lower-par[eInd], upper-par[eInd])
    par[eInd] <- par[eInd] + eps
    if(df==Inf) dens <- dmvnorm(eps, sigma=scale*Cmat, log=TRUE) - dmvnorm(0.0*eps, sigma=scale*Cmat, log=TRUE)
    else dens <- dmvt(eps, sigma=scale*Cmat, df=df, log=TRUE) - dmvt(0.0*eps, sigma=scale*Cmat, df=df, log=TRUE)
  } else if (fullPost & !is.null(object.sim$thetaSampList)) {
    if(is.null(thetaSamp)) thetaSamp <- length(object.sim$thetaSampList)
    parRow <- sample(1:nrow(object.sim$thetaSampList[[thetaSamp]]), 1, prob=object.sim$thetaSampList[[thetaSamp]][,1])
    par <- as.vector(object.sim$thetaSampList[[thetaSamp]][parRow,-c(1:3)])
    #print(parRow)
  } else par <- object.sim$par
  
  ###
  ### Process parameters for C++
  ###
  if (!is.null(err.mfX)) {
    theta.errX <- par[1:n.errX]
    Hmat <- exp(2 * err.mfX %*% theta.errX)
  } else Hmat <- rep(0.0, N)
  if (!is.null(err.mfY)) {
    theta.errY <- par[(n.errX + 1):(n.errX + n.errY)]
    Hmat <- cbind(Hmat,exp(2 * err.mfY %*% theta.errY))
  } else Hmat <- cbind(Hmat, Hmat)
  if(!is.null(rho)){
    Hmat = cbind(Hmat, sqrt(Hmat[,1])*sqrt(Hmat[,2])*rho)
  } else {Hmat = cbind(Hmat, rep(0,N))}
  Hmat[noObs==1,] = 0
  theta.mov <- par[(n.errX + n.errY + 1):(n.errX + n.errY + 2 * n.mov)]
  sig2 <- exp(2 * (mov.mf %*% theta.mov[1:n.mov]))
  b <- exp(mov.mf %*% theta.mov[(n.mov + 1):(2 * n.mov)])
  #activity <- rep(1, N)
  if (!is.null(activity)) {
    theta.activ <- par[(n.errX + n.errY + 2 * n.mov + 1)]
    b <- b / ((activity) ^ exp(theta.activ))
    active <- ifelse(b==Inf, 0, 1)
    b <- ifelse(b==Inf, 0, b) 
  } else active = rep(1,N)
  if (driftMod) {
    theta.drift <- par[(n.errX + n.errY + 2 * n.mov + 1):
                         (n.errX + n.errY + 2 * n.mov + 2)]
    b.drift <- exp(log(b) - log(1+exp(theta.drift[2])))
    sig2.drift <- exp(log(sig2) + 2 * theta.drift[1]) 
    out=CTCRWSAMPLE_DRIFT(y, Hmat, b, b.drift, sig2, sig2.drift, delta, noObs, active, a, P)
  } else {
    out=CTCRWSAMPLE(y, Hmat, b, sig2, delta, noObs, active, a, P)
  }
  
  if(driftMod){
    colnames(out$sim) <- apply(expand.grid(c("mu","theta","gamma"), c("x","y")), 1, paste, collapse=".")
  }  else {
    colnames(out$sim) <- apply(expand.grid(c("mu","nu"), c("x","y")), 1, paste, collapse=".")
  }
  ln.prior = ifelse(!is.null(object.sim$prior), object.sim$prior(par[eInd]), 0)
  isw <- ifelse(is.null(object.sim$thetaSampList) & fullPost==TRUE, out$ll - object.sim$loglik - dens, 0) + ln.prior
  samp <- list(alpha.sim=out$sim,
               locType=object.sim$locType, Time=object.sim$Time,
               loglik=out$lly+out$llx, par=par, log.isw = isw)
  class(samp) <- c("crwIS","list")
  attr(samp,"coord") <- object.sim$coord
  attr(samp,"random.drift") <- object.sim$driftMod
  attr(samp,"activity.model") <- !is.null(object.sim$activity)
  return(samp)
}

