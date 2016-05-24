makeViterbimat <-
function(obs, transitionmatrix, emissionmatrix, includeZeroState){
  numstates <- dim(transitionmatrix)[1]
  logv <- matrix(NA, nrow = length(obs), ncol = numstates)
  # Take logarithms of matrices to prevent underflow:
  logtrans <- log(transitionmatrix)
  logemi <- log(emissionmatrix)
  
  # Initial probabilities (t==0) of the states, copy number set equal to 1.
  logv[1,] <- log(0)
  if(includeZeroState){logv[1,2] <- log(1)}
  else{logv[1,1] <- log(1)}
  
  for (i in 2:length(obs)) 
  {
    for (l in 1:numstates)
    {
      # Find the probabilility, if we are in state l, of observing the number of reads in a window. We have to do obs[i]+1 to get from
      # the observation to the correct column (because 0 is included in possible observations.)
      statelprobcounti <- logemi[l,(obs[i]+1)]
      maxstate <- max(logv[(i-1),] + logtrans[l,])
      logv[i,l] <- statelprobcounti + maxstate
    }
  }
  return(logv)
}
