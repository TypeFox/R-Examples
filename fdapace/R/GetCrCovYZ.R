#' Functional Cross Covariance between longitudinal variable Y and scalar variable Z
#' 
#' Calculate the raw and the smoothed cross-covariance between functional
#' and scalar predictors using bandwidth bw or estimate that bw using GCV
#' 
#' @param Ly List of N vectors with amplitude information
#' @param Lt List of N vectors with timing information
#' @param Ymu Vector Q-1 Vector of length nObsGrid containing the mean function estimate (You can get that from FPCA)
#' @param bw Scalar bandwidth for smoothing the cross-covariance function (if NULL it will be automatically estimated)
#' @param Z Vector N-1 Vector of length N with the scalar function values
#' @param Zmu Scalar with the mean of Z (if NULL it will be automaticall estimated)
#' @param support Vector of unique and sorted values for the support of the smoothed cross-covariance function (if NULL it will be automatically estimated)
#' If the variables Ly1 is in matrix form the data are assumed dense and only the raw cross-covariance is returned.
#' @return A list containing:
#' \item{smoothedCC}{The smoothed cross-covariance as a vector}
#' \item{rawCC}{The raw cross-covariance as a vector }
#' \item{bw}{The bandwidth used for smoohting as a scaler}
#' \item{score}{The GCV score associated with the scalar used}
#' @examples
#' yList <- list( runif(5),  c(1:3), c(2:4), c(4))
#' tList <- list( c(1:5), c(1:3), c(1:3), 4)
#' Z = rep(4,4) # Constant vector so the covariance has to be zero.
#' sccObj = GetCrCovYZ(bw=1, Z= Z, Ly=yList, Lt=tList, Ymu=rep(4,5))
#' @references
#' \cite{Yang, Wenjing, Hans-Georg Mueller, and Ulrich Stadtmueller. "Functional singular component analysis." Journal of the Royal Statistical Society: Series B (Statistical Methodology) 73.3 (2011): 303-324}
#' @export

GetCrCovYZ <- function(bw = NULL, Z, Zmu = NULL, Ly, Lt = NULL, Ymu = NULL, support = NULL){
   
  # If only Ly and Z are available assume DENSE data
  if( is.matrix(Ly) && is.null(Lt) && is.null(Ymu) ){
    rawCC <- GetRawCrCovFuncScal(Ly = Ly, Z = Z)
    return ( list(smoothedCC = NULL, rawCC = rawCC, bw = bw, score = NULL) )  
  }
  # Otherwise assume you have SPARSE data
  if( is.null(Zmu) ){
    Zmu = mean(Z,na.rm = TRUE);
  }  
  # Get the Raw Cross-covariance 
  ulLt = unlist(Lt)
  if (is.null(support) ){
    obsGrid = sort(unique(ulLt))
  } else {
    obsGrid = support
  }
  # Check that the length of Z and the length of Ly are compatible
  if (length(Z) != length(Ly)){
    stop("Ly and Z are not compatible (possibly different number of subjects).")
  }
  
  if (is.null(Ymu)){
    stop("Ymu is missing without default.")
  }
  
  rawCC = GetRawCrCovFuncScal(Ly = Ly, Lt = Lt, Ymu = Ymu, Z = Z, Zmu = Zmu )

  # If the bandwidth is known already smooth the raw CrCov
  if( is.numeric(bw) ){
    smoothedCC <- smoothRCC(rawCC, bw, obsGrid )
    score = GCVgauss1D( smoothedY = smoothedCC, smoothedX = obsGrid, 
                        rawX = rawCC$tpairn, rawY = rawCC$rawCCov, bw = bw)
    return ( list(smoothedCC = smoothedCC, rawCC = rawCC, bw = bw, score = score) )
  # If the bandwidth is unknown use GCV to take find it
  } else {
    # Construct candidate bw's
    h0 = 2.0 * Minb( sort(ulLt), 2+1); # 2x the bandwidth needed for at least 3 points in a window
    r = diff(range(ulLt))    
    q = (r/(4*h0))^(1/9);   
    bwCandidates = sort(q^(0:19)*h0);
    # Find their associated GCV scores
    gcvScores = rep(Inf, length(bwCandidates))
    for (i in 1:length(bwCandidates)){
      smoothedCC <- try(silent=TRUE, smoothRCC(rawCC, bw = bwCandidates[i], xout = obsGrid ))
      if( is.numeric(smoothedCC) ){
        gcvScores[i] = GCVgauss1D( smoothedY = smoothedCC, smoothedX = obsGrid, 
                                   rawX = rawCC$tpairn, rawY = rawCC$rawCCov, bw = bwCandidates[i])
      }
    }
#browser()
    # Pick the one with the smallest score
    bInd = which(gcvScores == min(gcvScores, na.rm=TRUE));
    bOpt = max(bwCandidates[bInd]);
    smoothedCC <- smoothRCC( rawCC, bw = bOpt, obsGrid )
    return ( list(smoothedCC = smoothedCC, rawCC = rawCC, bw = bOpt, score = min(gcvScores, na.rm=TRUE)) )
  }  
}

# Calculate the smooth Covariances between functional and scalar predictors
# rCC     : raw cross covariance list object returned by GetRawCrCovFuncScal
# bw      : scalar
# xout    : vector M-1
# returns : vector M-1
smoothRCC <- function(rCC,bw,xout){
  x = matrix( unlist(rCC),  ncol=2)
  x= x[order(x[,1]),]
  return( Lwls1D(bw=bw, win=rep(1,nrow(x)), yin=x[,2], xin=x[,1], 'gauss', xout=xout) ) 
}

# Calculate GCV cost off smoothed sample assuming a Gaussian kernel
# smoothedY : vector M-1 
# smoothedX : vector M-1
# rawX      : vector N-1
# rawY      : vector N-1
# bw        : scalar
# returns   : scalar
GCVgauss1D <- function( smoothedY, smoothedX, rawX, rawY, bw){
  cvsum = sum( (rawY - approx(x=smoothedX, y=smoothedY, xout=rawX)$y)^2 );
  k0 = 0.398942;  # hard-coded constant for Gaussian kernel
  N  = length(rawX)
  r  = diff(range(rawX)) 
  return( cvsum / (1-(r*k0)/(N*bw))^2 )
}
