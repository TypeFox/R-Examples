stationdistTPM <-
function(M)  {
  eigenprob <- eigen(t(M))
  temp <- which(round(eigenprob$values,1)==1)
  stationdist <- eigenprob$vectors[,temp]
  stationdist <- stationdist/sum(stationdist)
  return(stationdist)
}
