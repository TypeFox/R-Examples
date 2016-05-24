
##==============================================================================
## plotellipse     : plots (part of) an ellipse
##==============================================================================

plotellipse <- function (rx=1, ry=0.2, mid=c(0,0), dr=0.01,
  angle=0, from=-pi, to=pi, type="l", lwd=2, lcol="black",
  col=NULL, arrow=FALSE, arr.length=0.4, arr.width=arr.length*0.5,
  arr.type="curved", arr.pos=1, arr.code=2, arr.adj=0.5,
  arr.col="black",  ...) {


  xy<-getellipse (rx,ry,mid,angle=angle,dr=dr,from=from,to=to)

  if (! is.null(col))
    polygon(xy,col=col,border=NA)
  if (type != "n" )
    lines(xy,type=type,lwd=lwd,col=lcol,...)
  nr <- nrow(xy)

  if (arrow) {
    ilen <- length(arr.pos)
    if (ilen>1) {
      arr.code  <- rep(arr.code  ,length.out=ilen)
      arr.col   <- rep(arr.col   ,length.out=ilen)
      arr.length<- rep(arr.length,length.out=ilen)
      arr.width <- rep(arr.width ,length.out=ilen)
      arr.type  <- rep(arr.type  ,length.out=ilen)
      arr.adj   <- rep(arr.adj   ,length.out=ilen)
    }

    for (i in 1: ilen) {
      ii <- max(2,trunc(nr*arr.pos[i]))
      Arrows(xy[ii-1,1], xy[ii-1,2], xy[ii,1], xy[ii,2],
            lcol=arr.col[i], code=arr.code[i], arr.col=arr.col[i],
            arr.length =arr.length[i], arr.width=arr.width[i],
            arr.type=arr.type[i], arr.adj=arr.adj[i])
    }
  }
}
