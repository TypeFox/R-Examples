
##==============================================================================
## filledcylinder  : draws and colors cylinder colors depend on radius
##==============================================================================

filledcylinder <- function(rx=1, ry=rx, len=1,
  col=femmecol(100), lcol=NA, lwd=2, lcolint=NULL,
  ltyint=1, lwdint=lwd, mid=c(0,0), angle=0, delt =1.0,
  dr=0.01, topcol=NULL, botcol=NULL, ...)  {

  rx       <- max(1e-6,rx)  # avoid NANs
  ncol     <- length (col)

  if (ncol > 0) {
    intrad   <- seq(pi/2,3*pi/2,length.out=ncol+1)
    Col      <- col
    nval <- ncol

    ## main body of cylinder
    for (i in 1:nval) {
      from <- intrad[i+1]
      to <- intrad  [i]

      cylindersegment (rx=rx,ry=ry,from=from,to=to, len=len,
                       mid=mid,dr=dr,angle=angle,
                       col=Col[i],delt=delt)
    }
  }

  base2 <- mid +c(len/2,0)
  if (angle != 0)
    base2 <- rotatexy(base2, angle=angle, mid=mid)

  base1 <- mid + c(-len/2,0)
  if (angle != 0)
    base1 <- rotatexy(base1, angle=angle, mid=mid)

  ## color of top
  if(! is.null (topcol))
    filledellipse( rx1=rx*delt, ry1=ry*delt, col=topcol, mid=base2,
                   angle=angle, dr=dr, ...)

  if (! is.null(botcol))
    filledellipse( rx1=rx, ry1=ry, col=botcol, mid=base1,
                   angle=angle, dr=dr, ...)


  if (! is.na(lcol)) {
    l1 <- rotatexy( getellipse( rx, ry, mid=mid+c(-len/2,0), dr=dr,
                    from=pi/2, to=3*pi/2), angle=angle, mid)
    l2 <- rotatexy( getellipse( rx*delt, ry*delt, mid=mid+c(len/2,0),dr=dr),
                    angle=angle, mid)

    lines(l1,col=lcol,lwd=lwd)
    lines(l2,col=lcol,lwd=lwd)
    if (! is.null(lcolint)) {
      l1 <- rotatexy(getellipse( rx, ry, mid=mid+c(-len/2,0), dr=dr,
                     from=-pi/2, to=pi/2), angle=angle, mid)
      lines(l1, col=lcolint, lwd=lwdint, lty=ltyint)
    }

    l1 <- rotatexy( rbind( mid+c(len/2, ry*delt), mid+c(-len/2,ry)),
                    angle, mid)
    l2 <- rotatexy( rbind( mid+c(len/2, -ry*delt), mid+c(-len/2,-ry)),
                    angle, mid)
    lines(l1,col=lcol,lwd=lwd)
    lines(l2,col=lcol,lwd=lwd)
  }

}

