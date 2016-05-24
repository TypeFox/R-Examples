
##==============================================================================
## getellipse   : calculates the x-y values for (part of) an ellipse
##==============================================================================

getellipse <- function (rx=1, ry=rx, mid=c(0,0), dr=0.01,
  angle=0, from=-pi, to=pi) {

  dr <- abs(dr)
  if (to < from) to <- 2*pi + to
  x  <- c( seq(from,to,by=dr), to)
  if (x[length(x)] == x[length(x)-1])
    x <- x[-length(x)]
  xy <- cbind( mid[1] + rx * cos(x), mid[2] + ry * sin(x))

  if (angle != 0)
    xy <- rotatexy (xy, angle=angle, mid=mid)  # rotate around mid
  return(xy)
}
