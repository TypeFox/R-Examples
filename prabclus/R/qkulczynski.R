qkulczynski <- function(regmat, log.distance=FALSE){
  if (log.distance)
    regmat <- log(regmat+1)
  nart <- ncol(regmat)
  jdist <- rep(0, nart*nart)
  dim(jdist) <- c(nart,nart)
  for (i in 1:(nart-1)){
#    cat("Row ",i,"\n")
    for (j in (i+1):nart){
      ri <- sum(regmat[,i])
      rj <- sum(regmat[,j])
      srij <- sum(pmin(regmat[,i],regmat[,j])) 
      jdist[j,i] <- jdist[i,j] <- 1 - 0.5* (srij/ri + srij/rj)
      if (is.na(jdist[i,j]))
        cat("qkulczynski warning! NA at i=",i,", j=", j,
            " srij=",srij," ri=",ri,"rj=",rj,"\n")
    }
  }
  jdist
}
