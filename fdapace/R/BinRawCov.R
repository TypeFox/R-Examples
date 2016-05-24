# Bin a `RawCov` object. Observations with the same time pairs are binned together and only one entry will be returned, containting the mean value (`meanVals`), weight (`count`), and residual sums of squares at each point (`RSS`). If `rcov$diag` is used then also bin the diagonal of the raw covariance similarly (with fields `diagMeans`, `diagRSS`, and `diagCount`.
# rcov: A `RawCov` object.
# returns: A list of class `BinnedRawCov`. 
BinRawCov <- function(rcov) {
  
  # Get the count, mean raw cov, and residual sum of squares at each pair of observed time points.
  tmp <- aggregate(rcov$cxxn, list(rcov$tPairs[, 1], rcov$tPairs[, 2]), 
    function(yy) c(RCPPmean(yy), length(yy), RCPPvar(yy) * (length(yy) - 1)))
  
  tPairs <- unname(as.matrix(tmp[, 1:2]))
  summaryDat <- tmp[, 3]
  meanVals <- summaryDat[, 1]
  count <- summaryDat[, 2]
  RSS <- summaryDat[, 3] # Residual sum of squares. For implementing GCV.
  RSS[is.na(RSS)] <- 0
  
  diagRSS <- diagCount <- diagMeans <- tDiag <- NULL
  if (!is.null(rcov$diag)) {
    tmp <- aggregate(rcov$diag[, 2], list(rcov$diag[, 1]), 
      function(yy) c(RCPPmean(yy), length(yy), RCPPvar(yy) * (length(yy) - 1)))
      
    tDiag <- tmp[, 1]
    diagSummary <- tmp[, 2]
    diagMeans <- diagSummary[, 1]
    diagCount <- diagSummary[, 2]
    diagRSS <- diagSummary[, 3]
    diagRSS[is.na(diagRSS)] <- 0
  }
  
  res <- list(tPairs=tPairs, meanVals=meanVals, RSS=RSS, 
              tDiag=tDiag, diagMeans=diagMeans, diagRSS=diagRSS, 
              count=count, 
              diagCount=diagCount, 
              error=rcov$error, dataType=rcov$dataType)
  class(res) <- 'BinnedRawCov'
  
  return(res)
}
