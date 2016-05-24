betaA <- function(n,eigenvaluesA,tPADmdemiY,DdemiPA,ddlmini,k,index0){
  prov <- rev(sumvalpr(k,n,rev(1-eigenvaluesA),n-index0+1,n-ddlmini+1))
  prov1 <- matrix(prov*as.vector(tPADmdemiY),n,1)
  return(DdemiPA%*%prov1)
}

sumvalpr <- function(k,n,valpr,index1,index0) {
  res <- rep(1,n)
  if (k>1) {
    n <- length(valpr)
    if (is.na(index1)&is.na(index0)) {
      res <- (1-valpr^k)/(1-valpr)
    } else {
      selecti <- NULL
      if (!is.na(index1)) selecti <- 1:index1
      if (!is.na(index0)) selecti <- c(selecti,index0:n)
      res[-selecti] <- (1-valpr[-selecti]^k)/(1-valpr[-selecti])
      if (!is.na(index1)) res[index1] <- k
      if (any(!is.finite(res[-selecti]))) {
        prov <- rep(1,length(valpr[-selecti]))
        kk <- 1
        while (kk <k) {
          prov <- prov+valpr[-selecti]^kk
          kk <- kk+1
        }
        res[-selecti] <- prov
      }
    }
  }
  return(res)
}
