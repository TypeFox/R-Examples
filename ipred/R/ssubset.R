
ssubset <- function(y, k, strat=TRUE) {
  if (!is.factor(y)) stop("y is not of class factor")
  N <- length(y)
  nlevel <- table(y)
  nindx <- list()
  indx <- 1:N
  outindx <- list()
  if (strat) {
    for (j in 1:length(nlevel))
      nindx <- c(nindx, list(indx[which(y == levels(y)[j])]))
    kmat <- kfoldcv(k, N, nlevel)
    for (i in 1:k) {
      sset <- kmat[,i]
      kindx <- c()
      for (j in 1:length(nlevel)) {
        if (i > 1)
          kindx <- c(kindx, nindx[[j]][(sum(kmat[j,
                     1:(i-1)])+1):sum(kmat[j,1:i])])
        else
          kindx <- c(kindx, nindx[[j]][1:kmat[j,1]])
      }
      kindx <- kindx[!is.na(kindx)]
      outindx <- c(outindx, list(kindx))
    }
    return(outindx)
  } else {
    kmat <- kfoldcv(k, N)
    nindx <- indx
    for (i in 1:k) { 
      if (i > 1)
        outindx <- c(outindx,
                  list(nindx[(sum(kmat[1:(i-1)])+1):sum(kmat[1:i])]))
      else
        outindx <- c(outindx, list(nindx[1:kmat[1]]))
    }
  }
  return(outindx)
}
