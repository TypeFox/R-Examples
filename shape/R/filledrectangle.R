
##==============================================================================
## filledrectangle : draws and colors rectangle color can be palette
##==============================================================================

filledrectangle <- function(mid=c(0,0), wx=1, wy=wx,
  col=femmecol(100), values=NULL, zlim=NULL,
  lwd=2, lcol=NA, angle = 0, ...) {

  vv     <- val2col(values,zlim,col)
  intrad <- vv$intrad
  Col    <- vv$Col

  ncol   <- length (Col)

  for (i in 1:ncol) {
    x1 <- mid[1] - wx/2
    x2 <- mid[1] + wx/2

    y1 <- mid[2] - wy/2  +wy*intrad[i]
    y2 <- mid[2] - wy/2  +wy*intrad[i+1]
    xy <- cbind(c(x1,x1,x2,x2),c(y1,y2,y2,y1))

    if (angle != 0)
      xy <- rotatexy (xy,angle=angle,mid=mid) # rotate around mid
    polygon(xy[,1],xy[,2],col=Col[i],border=Col[i], ...)

  }

  if (!is.na(lcol)) {
    x1 <- mid[1] - wx/2
    x2 <- mid[1] + wx/2

    y1 <- mid[2] - wy/2
    y2 <- mid[2] + wy/2
    xy <- cbind(c(x1,x1,x2,x2),c(y1,y2,y2,y1))

    if (angle != 0)
      xy <- rotatexy (xy,angle=angle,mid=mid)

    polygon(xy[,1],xy[,2],border=lcol,col=NA, ...)
  }
}
