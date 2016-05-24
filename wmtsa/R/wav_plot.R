################################################
## WMTSA package plot functionality
##
##  Functions:
##
##   wavStackPlot.default
##
################################################

###
# wavStackPlot.default
##

"wavStackPlot.default" <- function(x, x.axis=TRUE, y.axis=TRUE, type="l", plot=TRUE,
  bars=FALSE, vgap=.05, grid=FALSE, times=time(x[[1]]),
  grid.lty="dashed", same.scale=NULL,
  zerocenter=FALSE, zeroline=FALSE, col=rep(1,n), complex.math="mod", cex.main=0.7, cex.axis=0.7, ...)
{

  # develop local functions

  "ComplexMath" <- function(x, mode="mod")
  {
    if (is.complex(x))
      x <- switch(lowerCase(mode),
        mod=Mod(x), im =Im(x), arg=Arg(x), re=Re(x))
    x
  }

  if(is.numeric(x))
    x <- list(x)

  nms <- names(x)
  x   <- oldUnclass(x)
  n   <- length(x)

  if(is.null(nms))
    nms <- as.character(1:n)
  if(!is.element(length(type), c(1, n)))
    stop("type is the wrong length, should be", 1, "or", n)
  if(length(type) != n)
    type <- rep(type, n)
  ycenters <- mi <- rep(0., n)

  # compute yrange
  for(i in 1:n){
    yi <- ComplexMath(x[[i]], mode=complex.math)
    if(zerocenter)
      mi[[i]] <- 2*max(abs(yi), na.rm=TRUE)*(1+vgap)
    else
      mi[[i]] <- diff(range(yi, na.rm=TRUE))*(1+vgap)
    if(mi[[i]] < .Machine$double.eps*100)   # constant case
      mi[[i]] <- diff(range(c(0, yi), na.rm=TRUE))*(1+vgap)
  }

  if(length(same.scale)==1 && is.logical(same.scale)){
    if(same.scale)
      same.scale <- 1:n
    else
      same.scale <- numeric(0)
  }
  if(length(same.scale)){
    maxmi <- max(unlist(mi[same.scale]), na.rm=TRUE)
    for(i in same.scale)
      mi[[i]] <- maxmi
  }

  # compute xrange
  if (length(times) > 1)
    tt.delta <- times[2] - times[1]
  else
    tt.delta <- 0

  tt.start <- times[1] - tt.delta/2
  tt.range <- times[length(times)] - times[1] + tt.delta

  if (plot){
    xlim <- tt.start + c(0, tt.range) + c(-.03, .03+bars*.03)*tt.range
    ylim <- c(-n, 0)

    plot(xlim, ylim, type="n", axes=FALSE, xlab="", ylab="", xlim=xlim, ylim=ylim, cex.axis=cex.axis)
    if (x.axis)
      axis(side=1, at=pretty(times), line=1, srt=0, cex.axis=cex.axis, ...)
    if (y.axis)
      axis(side=2, at=-(1:n)+0.5, labels=nms, tick=FALSE, line=-1, srt=0, cex.axis=cex.axis, ...)
  }

  for(i in 1:n){
    yi <- ComplexMath(x[[i]], mode=complex.math)
    ni <- length(yi)
    tt <- ((0:(ni-1))/ni + 1/(2*ni))*tt.range + tt.start
    if(zerocenter){
      ycenters[[i]] <- -i+.5
      yi <- yi/mi[[i]]+ycenters[[i]]
      yzero <- TRUE
    }
    else if(all(yi>0, na.rm=TRUE) || all(yi<0, na.rm=TRUE)){
      ycenters[[i]] <- -i
      my <- min(yi, na.rm=TRUE)
      if(all(yi==my)) my <- 0
      yi <- (yi-my)/mi[[i]] + ycenters[[i]] + vgap/2
      yzero <- FALSE
    }
    else if(all(yi==0, na.rm=TRUE)){
      ycenters[[i]] <- -i+.5
      yi <- yi+ycenters[[i]]
      yzero <- TRUE
    }
    else{
      ycenters[[i]] <- -i + .5 - mean(range(yi, na.rm=TRUE)/mi[[i]])
      yi <- yi/mi[[i]] + ycenters[[i]]
      yzero <- TRUE
    }
    if(plot){
      if(type[i]=="h")
	segments(tt, rep(ycenters[[i]], ni), tt, as.vector(yi),col=col[i])
      else if(type[i]=="s"){
	tt <- tt - (tt[2]-tt[1])/2
	tt <- c(tt, tt[length(tt)]+(tt[2]-tt[1]))
	tt <- as.vector(rbind(tt, tt))
	yi <- as.vector(rbind(c(yi[2], yi), c(yi, yi[length(yi)-1])))
	lines(tt, yi, col=col[i], ...)
      }
      else lines(tt, yi, type=type[i], col=col[i], ...)
      if(yzero && zeroline)
	segments(tt[1], ycenters[[i]], tt[length(tt)], ycenters[[i]])
    }
  }
  if(plot && bars){
    mi <- 1/(mi)
    mi <- mi/max(mi)
    xloc <- times[length(times)] + .03*tt.range
    segments(xloc, -(1:n)+.5+mi/2, xloc, -(1:n)+.5-mi/2, lwd=2)
  }
  if(plot && grid) abline(h=-(1:(n-1)), lty=grid.lty)
#  xc <- sapply(x, length)
#  invisible(list(xc, ycenters))
  invisible(ycenters)
}
