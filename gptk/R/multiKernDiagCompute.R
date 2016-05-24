multiKernDiagCompute <-
function (kern, x) {
  if ( is.list(x) ) {
    xdim <- 0
    for ( i in seq_along(x) )
      xdim <- xdim + dim(as.array(x[[i]]))[1]

    k <- matrix(0, xdim, 1)
    startVal <- 1
    endVal <- dim(as.array(x[[1]]))[1]
    for ( i in seq(along=kern$comp) ) {

      k[startVal:endVal] <- kernDiagCompute(kern$comp[[i]], x[[i]])
      startVal <- endVal + 1
      if ( (i+1)<=length(kern$comp) )
        endVal <- endVal + dim(as.array(x[[i+1]]))[1]
    }
  } else {
    xdim <- dim(as.array(x))[1]
    k <- array(0, xdim*kern$numBlocks)
    startVal <- 1
    endVal <- xdim
    for ( i in seq(along=kern$comp) ) {
      k[startVal:endVal] <- kernDiagCompute(kern$comp[[i]], x)
      startVal <- endVal + 1
      endVal <- endVal + xdim
    }
  }

  return (k)
}
