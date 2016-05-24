
##==============================================================================
## plotcircle      : plots (part of) a circle
##==============================================================================

plotcircle <- function (r=1, ...) {

  user <- par("usr")
  pin  <- par("pin")
  sy   <- user[4] - user[3]
  sx   <- user[2] - user[1]
  ry   <- r * sy / sx * pin[1] / pin[2]

  plotellipse(rx=r, ry=ry, ... )
}
