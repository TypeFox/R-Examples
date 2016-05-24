multiKernExpandParam <-
function (kern, params) {
  if ( is.list(params) )
    params <- params$values

  params <- params%*%t(kern$paramGroups)
  startVal <- 1
  endVal <- 0

  for ( i in seq(along=kern$comp) ) {
    endVal <- endVal+kern$comp[[i]]$nParams
    kern$comp[[i]] <- kernExpandParam(kern$comp[[i]], params[startVal:endVal])
    startVal <- endVal+1
  }

  return (kern)
}
