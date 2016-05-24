
##==============================================================================
## filledellipse   : draws and colors ellipse; colors depend on radius
##==============================================================================

filledellipse <- function(rx1=1, rx2=0, ry1=rx1, ry2=NULL,
  mid=c(0,0), dr=0.01, angle=0, from=-pi, to=pi,
  col=femmecol(100), values=NULL, zlim=NULL, lwd=2, lcol=NA, ...) {

  ncol   <- length (col)
  if (is.null (ry2))
    ry2 <- ry1*rx2/rx1   # same proportionality as long radius

  ## outer ellipse
  xyouter   <- getellipse(rx=rx1, ry=ry1, mid=mid, dr=dr,
                          angle=angle, from=from, to=to)

  ## inner ellipse
  if (rx2==0)
    xyinner <- mid else {
    xyinner <- getellipse(rx=rx2, ry=ry2, mid=mid, dr=dr, angle=angle,
                          from=from, to=to)
    }
  filledshape(xyouter=xyouter, xyinner=xyinner, col=col, values=values,
              zlim=zlim, lwd=lwd, lcol=lcol, ...)

}
