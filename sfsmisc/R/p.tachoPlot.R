p.tachoPlot <- function(x, y, z, angle= c(pi/4,3*pi/4), size,
                        method= c("robust", "sensitive", "rank"),
                        legend = TRUE, show.method= legend,
                        xlab= deparse(substitute(x)),
                        ylab= deparse(substitute(y)), xlim, ylim, ...)
{
  ## Purpose: Puts a symbol (pointer) on a plot at each of the
  ##          specified locations.
  ## -------------------------------------------------------------------------
  ## Arguments: see on-line help (?p.tachoPlot)
  ## -------------------------------------------------------------------------
  ## Author: Christian Keller, Date: 16 Jun 95, 18:35

  if(length(angle) != 2)
    stop("length of angle must be 2")
  if(angle[1]<=0 | angle[1]>=pi/2)
    stop("angle[1] should be between 0 and pi/2")
  if(angle[2]<=pi/2 | angle[2]>=pi)
    stop("angle[2] should be between pi/2 and pi")

  method <- match.arg(method)

  xlab ; ylab ## eval substitute(.) now

  ii <- !is.na(x) & !is.na(y)
  x <- x[ii]; y <- y[ii]; z <- z[ii]

  if(method=="sensitive"){
    Min <- min(z, na.rm=TRUE)
    Max <- max(z, na.rm=TRUE)
    b <- (z-Min)/(Max-Min)
  }
  else if(method=="robust"){
    Range <- rrange(z)
    Min <- Range[1]
    Max <- Range[2]
    b <- pmin(pmax(z-Min,0),Max-Min)/(Max-Min)
  }
  else if(method=="rank"){
    Min <- min(z, na.rm=TRUE)
    Max <- max(z, na.rm=TRUE)
    Rank <- replace(rep(NA,length(z)), !is.na(z), rank(z[!is.na(z)]))
    b <- (Rank-1)/(sum(!is.na(z))-1)
  } else stop("unknown method (impossible)")

  ## -- range of the Plot
  range.x <- range(x)
  range.y <- range(y)
  pcm <- par("pin") * 2.54
  if(missing(size))
    size <- min(pcm)/20
  else {
    if(length(size) != 1)
      stop("length of size must be 1")
  }
  size <- size/2
  sx <- size*max(c(abs(cos(pi-angle[1])),abs(cos(pi-angle[2]))))
  sy <- size*max(c(abs(sin(pi-angle[1])),abs(sin(pi-angle[2]))))
  fx <- sx*diff(range.x)/(pcm[1]-2*size)
  fy <- sy*diff(range.y)/(pcm[2]-2*size)

  if(missing(xlim)) xlim <- range.x + c(-1,1)*fx
  if(missing(ylim)) ylim <- range.y + c(-1,1)*fy
  plot(x, y, pch=".", xlim=xlim, ylim=ylim, xlab=xlab, ylab=ylab, ...)

  ## -- calculate angles
  alpha <- angle[1] + (angle[2]-angle[1])*b
  usr <- par("usr")
  xd <- size*cos(pi-alpha)*diff(usr[1:2])/pcm[1]
  yd <- size*sin(pi-alpha)*diff(usr[3:4])/pcm[2]

  ## -- draw symbols
  if(method == "robust"){
    out <- z<Min | z>Max
    segments((x+xd)[!out],(y+yd)[!out], (x-xd)[!out], (y-yd)[!out], lty=1)
    if(any(out,na.rm=TRUE)) {
      segments((x+xd)[out],(y+yd)[out], (x-xd)[out], (y-yd)[out], lty=2,col=2)
    }
  }
  else{
    segments(x+xd, y+yd, x-xd, y-yd, lty=1)
  }
  if(legend){## -- draw legend
    cxy <- par("cxy")
    x1 <- min(pcm)/20*cos(pi-angle[1])*diff(usr[1:2])/pcm[1]
    x2 <- min(pcm)/20*cos(pi-angle[2])*diff(usr[1:2])/pcm[1]
    y1 <- min(pcm)/20*sin(pi-angle[1])*diff(usr[3:4])/pcm[2]
    y2 <- min(pcm)/20*sin(pi-angle[2])*diff(usr[3:4])/pcm[2]
    x <- usr[2] - 3*cxy[1] - x2
    y <- cxy[2] + usr[4]
    lines(c(x+x1,x,x+x2), c(y+y1,y,y+y2), lty=1, xpd=TRUE)
    text(x+x2, y, labels=formatC(Max), adj=0, cex=0.8*par("cex"), xpd=TRUE)
    text(x+x1, y, labels=formatC(Min), adj=1, cex=0.8*par("cex"), xpd=TRUE)
  }
  if(show.method)  ## -- print method name
    mtext(paste("method =",method),line=0, adj=1, cex=0.8*par("cex"))
  invisible()
}
