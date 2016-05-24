cmpndKernCompute <-
function (kern, x, x2) {
  if ( nargs()>2 ) {
    i <- 1
    if ( !is.na(kern$comp[[i]]$index) ) {
      k <- kernCompute(kern$comp[[i]], x[,kern$comp[[i]]$index], x2[,kern$comp[[i]]$index])
    } else {
      k <- kernCompute(kern$comp[[i]], x, x2)
    }
    for ( i in seq(2, length.out=(length(kern$comp)-1)) )
      if ( !is.na(kern$comp[[i]]$index) ) {
        k <- k+kernCompute(kern$comp[[i]], x[,kern$comp[[i]]$index], x2[,kern$comp[[i]]$index])
      } else {
        k <- k+kernCompute(kern$comp[[i]], x, x2)
      }
  } else {
    i <- 1
    if ( !is.na(kern$comp[[i]]$index) ) {
      k <- kernCompute(kern$comp[[i]], x[,kern$comp[[i]]$index])
    } else {
      k <- kernCompute(kern$comp[[i]], x)
    }
    for ( i in seq(2, length.out=(length(kern$comp)-1)) )
      if ( !is.na(kern$comp[[i]]$index) ) {
        k <- k+kernCompute(kern$comp[[i]], x[,kern$comp[[i]]$index])
      } else {
        k <- k+kernCompute(kern$comp[[i]], x)
      }
  }
  return (k)  
}
