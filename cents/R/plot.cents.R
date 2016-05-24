plot.cents <- function(x, y, xlab="t", ylab=expression(y[t]), marko=TRUE, ...){
  plot(x$y, type="l", xlab=xlab, ylab=ylab)
  T <- 1:length(x$y)
  Y <- x$censorPts[,1]
  segments(T, Y, T, Y, col=rgb(1,0,0,0.6))
  Y <- x$censorPts[,2]
  segments(T, Y, T, Y, col=rgb(1,0,0,0.6))
  if (marko) {
    xt <- 1:length(x$y)
    ind <- x$iy=="o"
    points(xt[ind],x$y[ind], pch=18, col=rgb(0,0,1,0.8))
    ind <- (x$iy=="L") | (x$iy=="R")
    points(xt[ind],x$y[ind], pch=1, col=rgb(1,0,1,0.8))
  }
}