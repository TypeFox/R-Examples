#####################
#' Robust covariance matrix estimation
#' @description Obtains a robust estimate of the covariance matrix of a sample of multivariate data, 
#'              using Campbell's (1980) method as described on p231-235 of Krzanowski (1988). 
#' @param sY A matrix, where each column is a replicate observation on a multivariate r.v.
#' @param alpha tuning parameter, see details.
#' @param beta tuning parameter, see details.
#' @details Campbell (1980) suggests an estimator of the covariance matrix which downweights observations 
#'                          at more than some Mahalanobis distance \code{d.0} from the mean.
#'                          \code{d.0} is \code{sqrt(nrow(sY))+alpha/sqrt(2)}. Weights are one for observations 
#'                          with Mahalanobis distance, \code{d}, less than \code{d.0}. Otherwise weights are 
#'                          \code{d.0*exp(-.5*(d-d.0)^2/beta)/d}. The defaults are as recommended by Campbell.
#'                          This routine also uses pre-conditioning to ensure good scaling and stable 
#'                          numerical calculations.
#' @return A list where:
#'         \itemize{
#'         \item{\code{E}}{a square root of the inverse covariance matrix. i.e. the inverse cov 
#'                         matrix is \code{t(E)\%*\%E};}
#'         \item{\code{half.ldet.V}}{Half the log of the determinant of the covariance matrix;}
#'         \item{\code{mY}}{The estimated mean;} 
#'         \item{\code{sd}}{The estimated standard deviations of each variable.}
#'          }
#' @references Krzanowski, W.J. (1988) Principles of Multivariate Analysis. Oxford.
#'             Campbell, N.A. (1980) Robust procedures in multivariate analysis I: robust covariance estimation. JRSSC 29, 231-237. 
#' @author Simon N. Wood, maintained by Matteo Fasiolo <matteo.fasiolo@@gmail.com>.

.robCov <- function(sY, alpha=2, beta=1.25, alpha2 = NULL, beta2 = NULL, tolVar = 10 * .Machine$double.eps) {
  ## Uses Campbell's robust approach as described on p 231 of Krzanowski 1988
  ## But adds pre-conditioning for stable computation....
  nStats <- nrow(sY)
  nObs <- ncol(sY)
  mY <- apply(sY, 1, median)
  sY1 <- sY - mY 
  ## use pre-conditioning to stabilize computation
  D <- rowMeans(sY1*sY1)^.5 
  
  # Identify statistics with very low variance and removing them
  lowVar <- which( D < tolVar )
  nLow <- length(lowVar)
  if( nLow )
  {
    nStats <- nStats - nLow
    if( nStats == 0 ) stop( paste("All the statistics have variance <", tolVar) )
    sY1 <- sY1[-lowVar, , drop = FALSE]
    mY <- mY[-lowVar]
    D <- D[-lowVar]
    warning( paste("There are", nLow, "statistics with variance <", tolVar, ": they were removed.") )
  }
  
  Di <- 1 / D  ## diagonal pre-conditioner
  
  sY1 <- Di*sY1 ## pre-conditioned for better scaling
  R <- qr.R(qr(t(sY1)))/sqrt(ncol(sY1)-1) ## Va = DR'RD - initial estimate
  zz <- forwardsolve(t(R),sY1)
  d <- sqrt(colSums((zz)^2)) ## Mahalonobis distance for each column
  
  ## create Campbell weight vector...
  d0 <- sqrt(nStats) + alpha/sqrt(2)
  w <- d*0 + 1
  ind <- d>d0
  w[ind] <- exp(-.5*(d[ind]-d0)^2/beta)*d0/d[ind] 
  
  ## The following lines all commented out because we want to use only the median for mY
  #if( nLow ) sY <- sY[-lowVar, , drop = FALSE]
  #mY <- colSums(w*t(sY))/sum(w)
  #sY1 <- sY - mY
  
  ## preconditioning...
  #D <- rowMeans(sY1*sY1)^.5
  #Di <- 1/D  ## diagonal pre-conditioner
  #sY1 <- Di*sY1 ## pre-conditioned for better scaling
  
  R <- qr.R(qr(w*t(sY1)))/sqrt(sum(w*w)-1) ## Va = DR'RD
  sd <- rowSums((D*t(R))^2)^.5
  E <- t(Di*backsolve(R,diag(nrow(R))))                   ## V^{-1} = E'E 
  half.ldet.V <- sum(log(abs(diag(R)))) + sum(log(D))
  COV <- tcrossprod(t(R)/Di, t(R)/Di) #Covariance matrix
  
  #Calculating weights of stats for stable saddlepoint
  if( !is.null(alpha2) && !is.null(beta2) )
  {
    d0 <- sqrt(nStats) + alpha2/sqrt(2)
    w <- d * 0 + 1
    ind <- d > d0
    w[ind] <- exp(-0.5 * (d[ind] - d0)^2/beta2) * d0/d[ind]
  }
  
  list(E=E, half.ldet.V=half.ldet.V, mY=mY, sd=sd, COV = COV, weights = w / sum(w) * nObs, lowVar = lowVar)
}
