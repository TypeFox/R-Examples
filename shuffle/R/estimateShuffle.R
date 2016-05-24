estimateShuffle <-
function(dat,prep,neg = FALSE){
  ########
  # Use the permutation variance estimator to estimate signal, noise and total variance
  # dat    - A vector/matrix of responses (each column a response)
  # prep   - A preprocessing structure, output of preparePVE
  # neg    - If true, allow negative variance estimates
  ########
  stopifnot(length(dat) == length(prep$des))

  # Get both original and permuted contrast
  Y_MSB  = MSbetAvg(dat, prep)
  PY_MSB = MSbetAvg(dat[prep$perm], prep)

  res = list()
  res$signalVar = prep$norm*(Y_MSB - PY_MSB)
  # Set negatives to -eps
  if (neg==FALSE & (res$signalVar < 0 )){
   res$signalVar = -0.00001
  }
  
  res$noiseVar =  Y_MSB - res$signalVar 
  res$YVar = Y_MSB
  res$effect =  res$signalVar/res$YVar
  
  return(res)
}
