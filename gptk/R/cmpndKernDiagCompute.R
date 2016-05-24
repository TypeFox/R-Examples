cmpndKernDiagCompute <-
function (kern, x) {
  i <- 1
  if ( !is.na(kern$comp[[i]]$index) ) {
    k <- kernDiagCompute(kern$comp[[i]], x[,kern$comp[[i]]$index])
  } else {
    k <- kernDiagCompute(kern$comp[[i]], x)
  }

  for ( i in seq(2, length=(length(kern$comp)-1)) )
    if ( !is.na(kern$comp[[i]]$index) ) {
      k <- k + kernDiagCompute(kern$comp[[i]], x[,kern$comp[[i]]$index])
    } else {
      k <- k + kernDiagCompute(kern$comp[[i]], x)
    }
       
  return (k)
}
