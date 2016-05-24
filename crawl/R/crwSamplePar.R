#' Create a weighted importance sample for posterior predictive track
#' simulation.
#' 

#' 
#' The \code{crwSamplePar} function uses a fitted model object from
#' \code{crwMLE} and a set of prediction times to construct a list from which
#' \code{\link{crwPostIS}} will draw a sample from either the posterior
#' distribution of the state vectors conditional on fitted parameters or a full
#' posterior draw from an importance sample of the parameters.
#' 

#' 
#' The crwSamplePar function uses the information in a
#' \code{\link{crwSimulator}} object to create a set of weights for importance
#' sample-resampling of parameters in a full posterior sample of parameters and
#' locations using \code{\link{crwPostIS}}. This function is usually called
#' from \code{\link{crwPostIS}}. The average user should have no need to call
#' this function directly.
#' 
#' @param object.sim A simulation object from \code{\link{crwSimulator}}.
#' @param method Method for obtaining weights for movement parameter samples
#' @param size Size of the parameter importance sample
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
#' @import mvtnorm
crwSamplePar <- function(object.sim, method="IS", size=1000, df=Inf, grid.eps=1, crit=2.5, scale=1)
{
  if(!inherits(object.sim, 'crwSimulator'))
    stop("Argument needs to be of class 'crwSimulator'\nUse 'crwSimulator( )' to create")
  fixPar <- object.sim$fixPar
  Cmat <- as.matrix(object.sim$Cmat[is.na(fixPar),is.na(fixPar)])
  se <- sqrt(diag(Cmat))
  err.mfX <- object.sim$err.mfX
  err.mfY <- object.sim$err.mfY
  parMLE <- object.sim$par
  n2ll.mode <- -2*object.sim$loglik
  activity <- object.sim$activity
  driftMod <- object.sim$driftMod
  mov.mf <- object.sim$mov.mf
  y <- object.sim$y
  noObs <- object.sim$noObs
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
  message("\nComputing importance weights ...\n")
  if(method=="IS"){
    thetaMat <- matrix(NA, size, length(fixPar)+3)
    for(i in 1:(size-1)){
      par <- parMLE
      eInd <- is.na(fixPar)
      eps <- rmvtt(mu=rep(0,sum(eInd)), Sigma=scale*Cmat, df=df, lower-par[eInd], upper-par[eInd])
      par[eInd] <- parMLE[eInd] + eps
      if(df==Inf) dens <- dmvnorm(eps, sigma=scale*Cmat, log=TRUE) - dmvnorm(0.0*eps, sigma=scale*Cmat, log=TRUE)
      else dens <- dmvt(eps, sigma=scale*Cmat, df=df, log=TRUE) - dmvt(0.0*eps, sigma=scale*Cmat, df=df, log=TRUE)
      ln.prior = ifelse(!is.null(prior), prior(par[eInd]), 0)
      n2ll.val <- crwN2ll(par[eInd], fixPar, y, noObs, delta, a,
                          P, mov.mf, err.mfX, err.mfY, rho=rho, activity=activity,
                          n.errX, n.errY, n.mov, driftMod, prior, need.hess=FALSE, 
                          constr=list(lower=lower, upper=upper)) + ln.prior
      thetaMat[i,] <- c(-n2ll.val/2 - dens, -n2ll.val/2, dens, par)
    }
    thetaMat[size,] <- c(object.sim$loglik, object.sim$loglik, 0, object.sim$par)
    thetaMat[,1] <- exp(thetaMat[,1]-max(thetaMat[,1]))/sum(exp(thetaMat[,1]-max(thetaMat[,1])))
  }
  else if(method=="mcmc"){
    
  }
  else if(method=="quadrature"){
    Eigen.list <- eigen(Cmat, symmetric=TRUE)
    V <- Eigen.list$vectors
    D <- diag(sqrt(Eigen.list$values))
    np <- sum(is.na(fixPar))
    grid.list <- rep(list(0), np)
    eInd <- is.na(fixPar)
    thetaMat <- matrix(c(-n2ll.mode/2, -n2ll.mode/2, 0, parMLE), nrow=1)
    for(k in 1:np){
      stop.grid <- TRUE
      z <- rep(0,np)
      while(stop.grid){
        z[k] <- z[k] + grid.eps 
        par <- parMLE
        par[eInd] <- parMLE[eInd] + V%*%D%*%z
        if(any(par[eInd]>upper) | any(par[eInd]<lower)) stop.grid <- FALSE
        else{
          n2ll.val <- crwN2ll(par[eInd], fixPar, y, noObs, delta, a,
                              P, mov.mf, err.mfX, err.mfY, rho=rho, activity=activity,
                              n.errX, n.errY, n.mov, driftMod, prior, need.hess=FALSE, 
                              constr=list(lower=lower, upper=upper))
          if(-(n2ll.mode - n2ll.val)/2 > crit) stop.grid <- FALSE
          else{
            grid.list[[k]] <- c(grid.list[[k]],z[k])
            thetaMat <- rbind(thetaMat, c(-n2ll.val/2, -n2ll.val/2, 0, par))
          }
        }
      }
      stop.grid <- TRUE
      z <- rep(0,np)
      while(stop.grid){
        z[k] <- z[k] - grid.eps 
        par <- parMLE
        par[eInd] <- parMLE[eInd] + V%*%D%*%z
        if(any(par[eInd]>upper) | any(par[eInd]<lower)) stop.grid <- FALSE
        else{
          n2ll.val <- crwN2ll(par[eInd], fixPar, y, noObs, delta, a,
                              P, mov.mf, err.mfX, err.mfY, rho=rho, activity=activity,
                              n.errX, n.errY, n.mov, driftMod, prior, need.hess=FALSE, 
                              constr=list(lower=lower, upper=upper))
          if(-(n2ll.mode - n2ll.val)/2 > crit) stop.grid <- FALSE
          else{
            grid.list[[k]] <- c(grid.list[[k]],z[k])
            thetaMat <- rbind(thetaMat, c(-n2ll.val/2, -n2ll.val/2, 0, par))
          }
        }
      }
    }
    grid.pts <- as.matrix(expand.grid(grid.list))
    grid.pts <- grid.pts[apply(grid.pts==0, 1, sum) < np-1, ]
    numEvals <- nrow(grid.pts)+nrow(thetaMat)
    message("\nEvaluating ", nrow(grid.pts)+nrow(thetaMat), " quadrature points ...\n")
    parFix <- ifelse(!eInd, parMLE, 0)
    for(i in 1:nrow(grid.pts)){
      z <- grid.pts[i,]
      par <- parMLE
      par[eInd] <- parMLE[eInd] + V%*%D%*%z
      if(any(par[eInd]>upper) | any(par[eInd]<lower)) next
      else{
        n2ll.val <- crwN2ll(par[eInd], fixPar, y, noObs, delta, a,
                            P, mov.mf, err.mfX, err.mfY, rho=rho, activity=activity,
                            n.errX, n.errY, n.mov, driftMod, prior, need.hess=FALSE, 
                            constr=list(lower=lower, upper=upper))
        if(-(n2ll.mode - n2ll.val)/2 > crit) next
        else thetaMat <- rbind(thetaMat, c(-n2ll.val/2, -n2ll.val/2, 0, par))
      }
    }
    thetaMat[,1] <- exp(thetaMat[,1]-max(thetaMat[,1]))/sum(exp(thetaMat[,1]-max(thetaMat[,1])))
  } else {
    stop("\nIncorrect specification of parameter sampling method\n")
  }
  
  colnames(thetaMat) <- c("w", "lik", "prop.lik", object.sim$nms)
  attr(thetaMat,"effSamp") <- nrow(thetaMat)/(1+(sd(thetaMat[,"w"])/mean(thetaMat[,"w"]))^2) 
  attr(thetaMat, "method") <- method
  attr(thetaMat, "numLikEval") <- ifelse(method=="quadrature", numEvals, size)
  if(is.null(object.sim$thetaSampList)) {
    object.sim$thetaSampList <- list(thetaMat)
  } else {
    object.sim$thetaSampList <- append(object.sim$thetaSampList, list(thetaMat))
  }
  return(object.sim)	
}

