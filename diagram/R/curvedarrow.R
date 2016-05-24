
##==============================================================================
# curvedarrow: Plot curved arrow at certain distance between two points
##==============================================================================

curvedarrow <- function(from, to, lwd=2, lty=1, lcol="black", arr.col=lcol, 
    arr.pos=0.5, curve=1, dr=0.01, endhead=FALSE, segment = c(0,1), ...)   {

  dpos  <- to-from
  angle <- atan(dpos[2]/dpos[1])*180/pi         # angle between both
  if (is.nan(angle)) return
  mid   <- 0.5*(to+from)                        # midpoint of ellipsoid arrow
  dst   <- dist(rbind(to, from))                # distance from-to
  ry    <- curve*dst                            # small radius of ellepsoid
  aFrom<-0                                      # angle to and from
  aTo<-pi
  if ( from[1] <= to[1]) {
    aFrom <- pi
    aTo <- 2*pi
  }

  if (segment [1] != 0)
    From <- segment[1] * aTo + (1-segment[1]) * aFrom
  else
    From <- aFrom
    
  if (segment [2] != 1)
    To <- segment[2] * aTo + (1-segment[2]) * aFrom
  else
    To <- aTo

  meanpi <- arr.pos * aTo + (1-arr.pos) * aFrom
  if (endhead) To <- meanpi


  plotellipse(rx=dst/2,  ry=ry, mid=mid, angle=angle, from = From, to = To,
             lwd=lwd, lty=lty, lcol=lcol)
  ell <- getellipse(rx=dst/2, ry=ry, mid=mid, angle=angle,
                    from=1.001*meanpi, to=0.999*meanpi, dr= 0.002)       #Changed from -0.002
  Arrows(ell[1,1], ell[1,2], ell[nrow(ell),1], ell[nrow(ell),2],
         code=1, lcol=lcol, arr.col=arr.col, ...)
  curvedarrow <- c(ell[nrow(ell),1], ell[nrow(ell),2])
}

