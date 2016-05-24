log.likelihood.calc <-
function(X, beta, death.times, ordered.time){
  log.likeli <- 0
  
  for(i in 1:length(death.times)){
    
    numer.ind <- which(ordered.time == death.times[i])
    denom.ind <- which(ordered.time >= death.times[i])
    numer <- sum(X[numer.ind,] %*% beta)
    denom <- (sum(exp(X[denom.ind,] %*% beta)))^length(numer.ind)

    log.likeli <- log.likeli + numer - log(denom)
  }
  return(-log.likeli)
}

