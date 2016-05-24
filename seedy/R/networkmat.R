networkmat <-
function(ID, sources) {
  tot.inf <- length(ID)
  truemat <- matrix(0, tot.inf, tot.inf)
  for (j in 1:tot.inf) {
    if (sources[j]%in%ID) {
      truemat[which(ID==sources[j]),j] <- 1
    }
  }
  return(truemat)
}
