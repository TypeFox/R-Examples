
##==============================================================================
# selfarrow: Plot circular arrow pointing at one point
##==============================================================================

selfarrow <- function(pos, lwd=2, lty=1, lcol="black", arr.pos=0.5,
    path="L", curve=c(0.1,0.1), dr=0.01, code=1, ...)   {

  if (length(curve)==1)
    curve <- c(curve, curve)
  mid <- pos
  if (path == "L") mid[1] <- mid[1] -curve[1]  # left shift
  if (path == "R") mid[1] <- mid[1] +curve[1]  # right shift
  if (path == "U") mid[2] <- mid[2] +curve[2]  # up shift
  if (path == "D") mid[2] <- mid[2] -curve[2]  # down shift

  aFrom <- 0
  aTo <- 2*pi                            # angle to and from

  meanpi <- arr.pos * aTo + (1-arr.pos) * aFrom       # pos of arrow
  plotellipse(rx=curve[1], ry=curve[2], mid=mid, from=aFrom, to=aTo,
             lwd=lwd, lty=lty, lcol=lcol)
  ell <- getellipse(rx=curve[1], ry=curve[2], mid=mid,
                    from=1.001*meanpi, to=0.999*meanpi, dr=-0.002)
  Arrows(ell[1,1], ell[1,2], ell[nrow(ell),1], ell[nrow(ell),2],
         code=code, lcol=lcol, ...)
  selfarrow <- c(ell[nrow(ell),1], ell[nrow(ell),2])
}
