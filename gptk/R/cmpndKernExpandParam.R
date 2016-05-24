cmpndKernExpandParam <-
function (kern, params) {
  if ( is.list(params) )
    params <- params$values
  params <- params %*% Matrix(t(kern$paramGroups))	## Params still log-transformed at this point.
  startVal <- 1
  endVal <- 0
  kern$whiteVariance <- 0
  for ( i in seq(along=kern$comp) ) {
    endVal <- endVal+kern$comp[[i]]$nParams
    kern$comp[[i]] <- kernExpandParam(kern$comp[[i]], params[seq(startVal,len=kern$comp[[i]]$nParams)])
    ## Used to be: kern$comp[[i]] <- kernExpandParam(kern$comp[[i]], params[startVal:endVal])
    startVal <- endVal+1
    if ( "white" %in% kern$comp[[i]]$type ) {
      kern$whiteVariance <- kern$whiteVariance+kern$comp[[i]]$variance
    } else if ( "whiteVariance" %in% names(kern$comp[[i]]) ) {
      kern$whiteVariance <- kern$whiteVariance+kern$comp[[i]]$whiteVariance
    }
  }      

  return (kern)
}
