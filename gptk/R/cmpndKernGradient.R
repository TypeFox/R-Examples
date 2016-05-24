cmpndKernGradient <-
function (kern, x, x2, covGrad) {
  if ( nargs()<4 ) 
    covGrad <- x2
  
  g <- array(0, dim(kern$paramGroups)[1])
  startVal <- 1
  endVal <- 0

  for ( i in seq(along=kern$comp) ) {
    endVal <- endVal + kern$comp[[i]]$nParams
    if ( !is.na(kern$comp[[i]]$index) ) {
      if ( nargs() < 4 ) {
	 g[seq(startVal,len=kern$comp[[i]]$nParams)] <- kernGradient(kern$comp[[i]], x[,kern$comp[[i]]$index], covGrad)
        ## Used to be: g[startVal:endVal] <- kernGradient(kern$comp[[i]], x[,kern$comp[[i]]$index], covGrad)
      } else {
        g[seq(startVal,len=kern$comp[[i]]$nParams)] <- kernGradient(kern$comp[[i]], x[,kern$comp[[i]]$index], x2[,kern$comp[[i]]$index], covGrad)
        ## Used to be: g[startVal:endVal] <- kernGradient(kern$comp[[i]], x[,kern$comp[[i]]$index], x2[,kern$comp[[i]]$index], covGrad)
      }
    } else {
      if ( nargs() < 4 ) {
        g[seq(startVal,len=kern$comp[[i]]$nParams)] <- kernGradient(kern$comp[[i]], x, covGrad)
        ## Used to be: g[startVal:endVal] <- kernGradient(kern$comp[[i]], x, covGrad)
      } else {
        g[seq(startVal,len=kern$comp[[i]]$nParams)] <- kernGradient(kern$comp[[i]], x, x2, covGrad)
        ## Used to be: g[startVal:endVal] <- kernGradient(kern$comp[[i]], x, x2, covGrad)
      }
    }
    startVal <- endVal + 1       
  }

  g <- g %*% kern$paramGroups    

  return (g)
}
