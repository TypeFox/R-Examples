qfsi <- function(nfactors, nstat, qscores, zsc_bn, qm) {
  #calculate FACTOR STABILITY INDEX
  mx <- sum(abs(qscores))*2 #maximum possible position changes
  #qm is the original analysis results
  fsi <- data.frame(FSindex=c(1:nfactors), NFSindex=c(1:nfactors))
  f <- 1
    while (f <= nfactors) {
        fsi[f,1] <- sum(abs(qm[[6]][f]-zsc_bn[f]))/nstat
        fsi[f,2] <- fsi[f,1]/mx
    f <- f+1
    }
  return(fsi)
}