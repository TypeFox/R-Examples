
##==============================================================================
## Arrows       : draws arrow with improved arrowhead
##==============================================================================

Arrows <- function(x0, y0, x1, y1, code=2,
  arr.length=0.4, arr.width=arr.length/2, arr.adj=0.5,
  arr.type="curved", segment=TRUE, col="black", lcol=col, lty=1,
  arr.col=lcol, lwd = 1, arr.lwd = lwd, ...)  {

  if (arr.type=="simple") {
    arrows(x0,y0,x1,y1,code=code,length=arr.length/2.54,
           lty=lty, col=col, lwd=lwd, ...)
    return()
  }
  if (arr.type=="T") {
    arrows(x0,y0,x1,y1,code=code,length=arr.length/(2*2.54),
           lty=lty, angle=90, col=col, lwd=lwd,  ...)
    return()
  }

  ## draw segment
  if (segment)                                # version 1.4: added lwd 
    segments(x0,y0,x1,y1,col=lcol,lty=lty,lwd=lwd,...)

  ## scaling factor
  user<-par("usr")
  pin <-par("pin")
  pin <- pin/max(pin)
  sy<- (user[4]-user[3]) /pin[2]
  sx<- (user[2]-user[1]) /pin[1]

  ## code = 2
  angle<- atan((y1-y0) /(x1-x0) *sx/sy)/pi*180
  angle[is.nan(angle)]<-0
  angle [x1<x0] <-180+angle[x1<x0]
  xx<-x1
  yy<-y1

  ## code =3 draws two arrowheads
  if (code == 3)
    Arrowhead(x0=xx,y0=yy,angle=angle,
              lcol=lcol,arr.col=arr.col,arr.adj=arr.adj,
              lty=lty,arr.length=arr.length,arr.width=arr.width,
              arr.type=arr.type,arr.lwd=arr.lwd)

  if (code != 2) {
    angle <-180 + angle
    xx<-x0
    yy<-y0
  }

  Arrowhead(x0=xx,y0=yy,angle=angle,lcol=lcol,arr.col=arr.col,
            arr.adj=arr.adj,lty=lty,arr.length=arr.length,
            arr.width=arr.width,arr.type=arr.type,arr.lwd=arr.lwd)
}


