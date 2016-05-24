expTransform <-
function (x, transform="atox") {

  eps <- 2.2204e-16
  maxVal <- 4.3112e+15

  thre <- 36	## threshold
  y <- array(0, dim(as.array(x)))

  if ( "atox" == transform ) {
    for ( ind in seq_along(as.array(x)) ) {
      if ( x[ind] > thre ) y[ind] <- maxVal else
      if ( x[ind] < -thre ) y[ind]<- eps else
      y[ind] <- exp(x[ind])
    }
  } else if ( "xtoa" == transform ) {
    for ( ind in seq_along(as.array(x)) ) {
      y[ind] <- .complexLog(x[ind])
    }
  } else if ( "gradfact" == transform )
    y <- x

  return (y)
}
