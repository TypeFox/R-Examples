
##==============================================================================
## filledmultigonal : draws and colors shape with equal-sized vertices
## color can be palette
##==============================================================================

filledmultigonal <- function(mid=c(0,0), rx=1, ry=rx, nr=4,
  col=femmecol(100), values=NULL, zlim=NULL, lwd=2, lcol=NA,
  angle = 0, ...) {

  xy <- getellipse(mid=mid,rx=rx,ry=ry,dr=2*pi/nr,from=0,to=2*pi)
  if (angle != 0)
    xy <- rotatexy (xy,angle=angle,mid=mid) # rotate around mid

  filledshape(xyouter=xy, xyinner=mid, col=col, values=values,
              zlim=zlim, lwd=lwd, lcol=lcol, ...)

}
