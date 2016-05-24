
##==============================================================================
##  ROTATIONS...
## rotatexy        : rotates xy values around midpoint
##==============================================================================

rotatexy   <- function (xy, angle, mid=colMeans(xy), asp=FALSE) {

  xy    <- matrix(ncol=2,data=xy)
  angpi <- angle / 180 *pi
  cosa  <-cos(angpi)
  sina  <-sin(angpi)

  dx    <- xy[,1] - mid[1]
  dy    <- xy[,2] - mid[2]

  ex    <-mid[1] + cosa*dx-sina*dy
  ey    <-mid[2] + sina*dx+cosa*dy

  if (asp) {
    user <- par("usr")
    pin  <- par("pin")
    sy   <- user[4]-user[3]
    sx   <- user[2]-user[1]
    ey   <- mid[2] + (ey -mid[2])*sy/sx*pin[1]/pin[2]
  }

  return(cbind(ex,ey))
}
