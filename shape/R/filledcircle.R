
##==============================================================================
## filledcircle    : draws and colors circle; colors depend on radius
##==============================================================================

filledcircle <- function(r1=1, r2=0, mid=c(0,0), dr=0.01,
  from=-pi, to=pi, col=femmecol(100), values=NULL, zlim=NULL,
  lwd=2, lcol=NA, ...) {

  user <- par("usr")             # to maintain the aspect ratio...
  pin  <- par("pin")
  sy   <- user[4]-user[3]
  sx   <- user[2]-user[1]
  ry1  <- r1*sy/sx*pin[1]/pin[2]
  ry2  <- r2*sy/sx*pin[1]/pin[2]
  filledellipse(rx1=r1,ry1=ry1,rx2=r2,ry2=ry2,mid=mid,dr=dr,from=from,to=to,
               col=col,lwd=lwd,lcol=lcol,values=values,zlim=zlim,...)
}
