#############################################################
#
#	totalvariation.circular
#       GNU General Public Licence 2.0 
#	Author: Claudio Agostinelli
#	E-mail: claudio@unive.it
#	Date: September, 17, 2012
#	Version: 0.6
#
#	Copyright (C) 2012 Claudio Agostinelli
#
#############################################################

totalvariation.circular <- function(x, y, z=NULL, q=0.95, bw, adjust = 1, type = c("K", "L"), kernel = c("vonmises", "wrappednormal"), na.rm = FALSE, step=0.001, eps.lower=10^(-4), eps.upper=10^(-4), ...) {
  if (is.null(z))
    z <- circular(seq(0,2*pi+step,step))
  x <- conversion.circular(x, units="radians", zero=0, rotation="counter", modulo="2pi")
  y <- conversion.circular(y, units="radians", zero=0, rotation="counter", modulo="2pi")
  z <- conversion.circular(z, units="radians", zero=0, rotation="counter", modulo="asis")  
  modalx <- modal.region.circular(x=x, z=z, q=q, bw=bw, adjust=adjust, type=type, kernel=kernel, na.rm=na.rm, step=step, eps.lower=eps.lower, eps.upper=eps.upper, ...)
  modaly <- modal.region.circular(x=y, z=z, q=q, bw=bw, adjust=adjust, type=type, kernel=kernel, na.rm=na.rm, step=step, eps.lower=eps.lower, eps.upper=eps.upper, ...)
  zerosx <- modalx$zeros
  areasx <- modalx$areas$tot
  zerosy <- modaly$zeros
  nx <- nrow(zerosx)
  ny <- nrow(zerosy)
  areasy <- modaly$areas$tot
  denx <- modalx$density
  deny <- modaly$density  
  denx <- denx$y/areasx
  deny <- deny$y/areasy
  denx[z < zerosx[1,1]] <- 0
  if (nx > 1) {
    for (i in 1:(nx-1)) {
      denx[z > zerosx[i,2] & z < zerosx[i+1,1]] <- 0
    }
  }
  denx[z > zerosx[nx, 2]] <- 0
  deny[z < zerosy[1,1]] <- 0
  if (ny > 1) {
    for (i in 1:(ny-1)) {
      deny[z > zerosy[i,2] & z < zerosy[i+1,1]] <- 0
    }
  }
  deny[z > zerosy[ny, 2]] <- 0  
  den <- ifelse(denx > deny, denx-deny, 0) 
  denmax <- approxfun(x=z, y=den)
  tv <- 0
##  byhand <- 0    
  for (i in 1:nx) {
##    byhand <- byhand + step*sum(den[z >= zerosx[i,1] & z <= zerosx[i,2]])
    tv <- tv +  integrate(denmax, lower=zerosx[i,1]+step, upper=zerosx[i,2]-step)$value
  }
  result <- list()
  result$tv <- tv
  result$ovl <- 1 - tv
##  result$byhand <- byhand
  result$q <- q
  result$bw <- bw
  result$modal.x <- modalx
  result$modal.y <- modaly
  result$density.x <- approxfun(x=z, y=denx)
  result$density.y <- approxfun(x=z, y=deny)
  result$density <- denmax
  class(result) <- 'totalvariation.circular'
  return(result)
}

plot.totalvariation.circular <- function(x, tv=TRUE, ovl=TRUE, units=c('radians', 'degrees', 'hours'), xlab=NULL, ylab=NULL, main=NULL, from=0, to=2*pi, add=FALSE, n=1000, polygon.control.tv=list(), polygon.control.ovl=list(), xlty=1, ylty=1, xcol=1, ycol=1, xlwd=1, ylwd=1, axes=TRUE, ...) {
  units <- match.arg(units)
  if (is.null(xlab))
    xlab <- paste('bw=', round(x$bw, 3), sep='')
  if (is.null(ylab))
    ylab <- 'Kernel density estimates'
  if (is.null(main))
    main <- 'Common area under the curves'  
  polygon.control.default <- list(density = NULL, angle = 45, border = NA, col = NA, lty = par("lty"), fillOddEven = FALSE)
  npc.tv <- names(polygon.control.tv)
  npc.ovl <- names(polygon.control.ovl)  
  npcd <- names(polygon.control.default)
  polygon.control.tv <- c(polygon.control.tv, polygon.control.default[setdiff(npcd, npc.tv)])
  polygon.control.ovl <- c(polygon.control.ovl, polygon.control.default[setdiff(npcd, npc.ovl)])  
  if (add==TRUE)
    plot(x$density.x, from=from, to=to, add=TRUE, n=n, xlab=xlab, ylab=ylab, main=main, lty=xlty, col=xcol, lwd=xlwd, ...)
  else {
    plot(x$density.x, from=from, to=to, add=FALSE, axes=FALSE, n=n, xlab=xlab, ylab=ylab, main=main, lty=xlty, col=xcol, lwd=xlwd, ...)
    if (axes) {
      axis(2)
      if (units=='degrees') {
        labels <- c(0, 45, 90, 135, 180, 225, 270, 315, 360)
        at <- labels*pi/180
      } else if (units=='hours') {
        labels <- c(0, 3, 6, 9, 12, 15, 18, 21, 24)
        at <- labels*pi/12
      } else {
        labels <- NULL
        at <- axTicks(1)
      }
      axis(1, at=at, labels=labels)
    }
  }
  plot(x$density.y, from=from, to=to, add=TRUE, n=n, lty=ylty, col=ycol, lwd=ylwd, ...)
  z <- seq(from, to, length.out=n)
  denmin <- pmin(x$density.x(z), x$density.y(z))
  denmax <- pmax(x$density.x(z), x$density.y(z))  
  if (tv)
    polygon(x=c(from, z, to, rev(z), from), y=c(0, denmin, 0, rev(denmax), 0), density = polygon.control.tv$density, angle = polygon.control.tv$angle, border = polygon.control.tv$border, col = polygon.control.tv$col, lty = polygon.control.tv$lty, fillOddEven = polygon.control.tv$fillOddEven)
  if (ovl)
    polygon(x=c(from, z, to, from), y=c(0,denmin, 0, 0), density = polygon.control.ovl$density, angle = polygon.control.ovl$angle, border = polygon.control.ovl$border, col = polygon.control.ovl$col, lty = polygon.control.ovl$lty, fillOddEven = polygon.control.ovl$fillOddEven)
  abline(v=c(x$modal.x$zeros), lty=2)
  abline(v=c(x$modal.y$zeros), lty=2)  
  invisible(x)
}

