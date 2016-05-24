t2p <- function (T) {
  if (is.null(dim(T))) {
    T <- array(T, dim=c(length(T), 1))
  }
  otherDims <- 2:length(dim(T))
  pvalue <- function (Tvec) {
    t0 <- Tvec[1] ## observed
    Tdist <- Tvec[-1] ## distribution
    B <- sum(!is.na(Tdist)) ## number of (valid) permutations
    r0 <- sum(Tdist >= t0, na.rm=TRUE)
    Rdist <- B - rank(Tdist, na.last='keep', ties.method='min') + 1
    Pvec <- c(r0, Rdist) / B
    return(Pvec) ## may include NAs
  }
  P <- apply(T, otherDims, pvalue)
  return(P)
}
