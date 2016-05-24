##==============================================================================
## cylindersegment      : plots (part of) a cylinder
##==============================================================================


cylindersegment <- function (rx=1, ry=rx, from=pi, to=3*pi/2,
  len=1, mid=c(0,0), angle=0, dr=0.01, col="black", delt =1.0, ...) {

  base1 <- mid ; base1[1]<-mid[1]-len/2
  base2 <- mid ; base2[1]<-mid[1]+len/2

  if (from >to)  {
    a <- from
    from <- to
    to <- a
  }

  x <- c( seq(from,to,by=dr), to)
  ex <- base1[1] + rx * cos(x)
  ey <- base1[2] + ry * sin(x)

  if (angle != 0) {
    xy <- rotatexy (cbind (ex,ey) , angle=angle, mid=mid)
    ex <- xy[,1];ey <- xy[,2]
  }

  x<-rev(x)
  x2 <- base2[1] + rx * cos(x) * delt
  y2 <- base2[2] + ry * sin(x) * delt

  if (angle != 0) {
    xy <- rotatexy ( cbind(x2,y2), angle=angle, mid=mid)
    x2 <- xy[,1];y2 <- xy[,2]
  }

  ex <- c(ex, x2)
  ey <- c(ey, y2)

  polygon(ex, ey, col=col, border=col, ...)

} # end cylindersegment
