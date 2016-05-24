holm <- function(pvals){
  p <- dim(pvals)[1]
  pvalsH <- pvals
  temp <- pvals[lower.tri(pvals)]
  oo <- order(temp)
  temp <- temp[oo]
  q <- p*(p-1)/2
  tempHa <- 1-(1-temp)^((q:1)/q)
  tempH <- rep(0,q)
  for(i in 1:q){
    tempH[i] <- max(tempHa[1:i])
  }
  tempH[oo] <-  tempH
  pvalsH[lower.tri(pvalsH)] <- tempH
  pvalsH[upper.tri(pvalsH)] <- t(pvalsH)[upper.tri(pvalsH)]
  return(zapsmall(pvalsH))
}
