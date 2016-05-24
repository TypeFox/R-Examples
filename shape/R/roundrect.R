
##==============================================================================
## roundrect: rectangular box with rounded left and right edges
##==============================================================================

roundrect<- function(mid, radx, rady, rx=rady, dr=0.01, col="white",
  lcol="black", lwd=2, angle=0, ...) {

  leftell  <- getellipse(ry=rady, rx=rx, mid=c(mid[1]-radx,mid[2]),
                         dr=dr, from=pi/2, to=3/2*pi)
  rightell <- getellipse(ry=rady, rx=rx, mid=c(mid[1]+radx,mid[2]),
                         dr=dr, from=-pi/2, to=pi/2)
  xyout <- rbind(leftell, c(mid[1]-radx,mid[2]-rady),
                 c(mid[1]+radx, mid[2]-rady),
                 rightell, c(mid[1]+radx,mid[2]+rady),
                 c(mid[1]-radx, mid[2]+rady) )
  if (angle != 0)
    xyout <- rotatexy(xyout, angle=angle)
  filledshape( xyout, mid, col=col, lcol=lcol, lwd=lwd, ...)

}
