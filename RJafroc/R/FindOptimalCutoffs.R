FindOptimalCutoffs <- function(nl, ll) {
  UNINITIALIZED <- RJafrocEnv$UNINITIALIZED
  nl <- nl[nl != UNINITIALIZED]
  nlLen <- length(nl)
  ll <- ll[ll != UNINITIALIZED]
  llLen <- length(ll)
  cutoff1 <- min(c(nl, ll)) - 0.1
  cutoffMax <- max(c(nl, ll)) + 0.1
  minBinCount <- 5
  maxCutoffs <- 15
  iniCutoff <- cutoff1
  
  cutoffFound <- FALSE
  while (!cutoffFound) {
    cutoff1 <- iniCutoff
    cutoff2 <- iniCutoff
    zeta <- iniCutoff
    orgCtffDelta <- (cutoffMax - cutoff1)/1000
    ctffDelta <- orgCtffDelta
    while (1) {
      if (cutoff1 >= cutoffMax) {
        cutoffFound <- TRUE
        break
      }
      nlCount <- sum(nl > cutoff2 & nl <= cutoff1)
      llCount <- sum(ll > cutoff2 & ll <= cutoff1)
      if (nlCount >= minBinCount && llCount >= minBinCount) {
        if (nlCount > minBinCount && llCount > minBinCount) {
          if (ctffDelta > 0) 
            ctffDelta <- -ctffDelta/2
          cutoff1 <- cutoff1 + ctffDelta
          if (abs(ctffDelta) < abs(orgCtffDelta)/256) {
            cutoff1 <- cutoff1 - ctffDelta
            zeta <- c(zeta, cutoff1)
            ctffDelta <- orgCtffDelta
            cutoff2 <- cutoff1
          }
        } else {
          zeta <- c(zeta, cutoff1)
          ctffDelta <- orgCtffDelta
          cutoff2 <- cutoff1
        }
        if (length(zeta) == maxCutoffs) {
          minBinCount <- minBinCount + 1
          break
        }
      } else {
        if (ctffDelta < 0) 
          ctffDelta <- -ctffDelta/2
        cutoff1 <- cutoff1 + ctffDelta
      }
    }
  }
  if (length(zeta) > 1) 
    zeta <- zeta[-length(zeta)]
  return(zeta)
} 
