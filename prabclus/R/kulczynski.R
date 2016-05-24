"kulczynski" <-
function(regmat){
  nart <- ncol(regmat)
  jdist <- rep(0, nart*nart)
  dim(jdist) <- c(nart,nart)
  for (i in 1:(nart-1)){
#    cat("Row ",i,"\n")
    for (j in (i+1):nart){
      ri <- sum(regmat[,i])
      rj <- sum(regmat[,j])
      srij <- sum(regmat[,i]+regmat[,j]>=1.99) 
      jdist[j,i] <- jdist[i,j] <- 1 - 0.5* (srij/ri + srij/rj)
      if (is.na(jdist[i,j]))
        cat("Warning! NA at i=",i,", j=", j,"\n")
    }
  }
  jdist
}
