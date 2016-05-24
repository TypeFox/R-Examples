trimmedMean <- function(pos, score, probeWindow=600, meanTrim=.1, nProbes=10) {
  st <- 1
  en <- 1
  n <- length(pos)
  stopifnot( length(score)==n )
  tmean <- rep(0,n)
  for(ii in 1:n) {
    while( (pos[ii]-pos[st]) > probeWindow )
      st <- st + 1
    while( (pos[en]-pos[ii]) < probeWindow & (en < n))
      en <- en + 1
    if ( (en-st+1) < nProbes )
      next
    lo <- floor((en-st+1) * meanTrim) + 1
    hi <- (en-st+1) + 1 - lo
	#cat(lo,hi,(en-st+1),"\n")
    tmean[ii] <- mean( score[st:en], trim=meanTrim )*sqrt(hi-lo+1)
  }
  tmean
}
