modal.region <- function(x, ...) UseMethod("modal.region")

modal.region.default <- function(x, ...) .NotYetImplemented()

#############################################################
#
#	modal.region.circular
#       GNU General Public Licence 2.0 
#	Author: Claudio Agostinelli
#	E-mail: claudio@unive.it
#	Date: July, 21, 2011
#	Version: 0.6
#
#	Copyright (C) 2011 Claudio Agostinelli
#
#############################################################

modal.region.circular <- function(x, z=NULL, q=0.95, bw, adjust = 1, type = c("K", "L"), kernel = c("vonmises", "wrappednormal"), na.rm = FALSE, step=0.01, eps.lower=10^(-4), eps.upper=10^(-4), ...) {
  if (is.null(z))
    z <- circular(seq(0,2*pi+step,step))
  if (is.circular(x))
    xcp <- circularp(x)     
  else
    xcp <- list(type="angles", units="radians", template="none", modulo="asis", zero=0, rotation="counter")   
  x <- conversion.circular(x, units="radians", zero=0, rotation="counter", modulo="2pi")
  z <- conversion.circular(z, units="radians", zero=0, rotation="counter", modulo="asis")  
  object <- density.circular(x=x, z=z, bw=bw, adjust=adjust, type=type, kernel=kernel, na.rm=na.rm)
  
  if (q > 1 | q < 0)
    stop("'q' must be between 0 and 1")
  pos.max <- which.max(object$y)
  pos.min <- which.min(object$y)  
  if (q==0 & length(pos.max)==1) {
    zeros <- matrix(c(object$x[pos.max], object$x[pos.max]), nrow=1)
    areas <- list(tot=0, areas=0)
    l <- object$y[pos.max]
  } else if (q==1 & length(pos.min)==1) {
    zeros <- matrix(c(0, object$x[pos.min], object$x[pos.min], 2*pi), nrow=2, byrow=TRUE)
    areas <- list(tot=1, areas=1)
    l <- object$y[pos.min]
  } else {
    internal <- function(x) {
      area(allcrosses(l=x, object=object, grid=z), object)$tot - q
    }
    l <- uniroot(internal, lower=min(object$y)*(1+eps.lower), upper=max(object$y)*(1-eps.upper))$root
    zeros <- allcrosses(l=l, object=object, grid=z)
    areas <- area(zeros, object)
  }
  result <- list()
  xunits <- circularp(x)$units
  result$zeros <- conversion.circular(circular(zeros), xcp$units, xcp$type, xcp$template, 'asis', xcp$zero, xcp$rotation)
  result$areas <- areas
  object$x <- conversion.circular(object$x, xcp$units, xcp$type, xcp$template, 'asis', xcp$zero, xcp$rotation)
  result$density <- object
  result$q <- q
  result$level <- l
  class(result) <- 'modal.region.circular'
  return(result)
}

#############################################################
#
#	plot.modal.region.circular
#       GNU General Public Licence 2.0 
#	Author: Claudio Agostinelli
#	E-mail: claudio@unive.it
#	Date: July, 5, 2011
#	Version: 0.5-1
#
#	Copyright (C) 2011 Claudio Agostinelli
#
#############################################################

plot.modal.region.circular <- function(x, plot.type=c('line', 'circle'), xlab=NULL, ylab=NULL, xlim=NULL, ylim=NULL, main=NULL, polygon.control=list(), ...) {
  polygon.control.default <- list(density = NULL, angle = 45, border = NULL, col = NA, lty = par("lty"), fillOddEven = FALSE)
  npc <- names(polygon.control)
  npcd <- names(polygon.control.default)
  polygon.control <- c(polygon.control, polygon.control.default[setdiff(npcd, npc)])
  plot.type <- match.arg(plot.type)
  if (is.null(xlab))
    xlab <- paste('bw=', round(x$density$bw,3), sep='')
  if (is.null(ylab))
    ylab <- 'Kernel Density Estimates'
  if (is.null(main))
    main <- 'Areas under the curve'
  plot(x$density, plot.type=plot.type, xlab=xlab, ylab=ylab, main=main, xlim=xlim, ylim=ylim, ...)
  if (plot.type=='line') {
    abline(h=x$level, lty=2)
    abline(v=c(x$zeros), lty=2)
    for (i in 1:nrow(x$zeros)) {
      zero1 <- x$zeros[i,1]
      zero2 <- x$zeros[i,2]
      inside <- x$density$x >= zero1 & x$density$x <= zero2
      polygon(x=c(zero2, zero1, x$density$x[inside], zero2), y=c(0,0,x$density$y[inside], 0), density = polygon.control$density, angle = polygon.control$angle, border = polygon.control$border, col = polygon.control$col, lty = polygon.control$lty, fillOddEven = polygon.control$fillOddEven)
    }
  } else {
    warning('Not Yet Implemented for plot.type=circle')
  }
}

#############################################################
#
#	lines.modal.region.circular
#       GNU General Public Licence 2.0 
#	Author: Claudio Agostinelli
#	E-mail: claudio@unive.it
#	Date: July, 5, 2011
#	Version: 0.5-1
#
#	Copyright (C) 2011 Claudio Agostinelli
#
#############################################################

lines.modal.region.circular <- function(x, plot.type=c('line', 'circle'), polygon.control=list(), ...) {
  polygon.control.default <- list(density = NULL, angle = 45, border = NULL, col = NA, lty = par("lty"), fillOddEven = FALSE)
  npc <- names(polygon.control)
  npcd <- names(polygon.control.default)
  polygon.control <- c(polygon.control, polygon.control.default[setdiff(npcd, npc)])  
  plot.type <- match.arg(plot.type)
  lines(x$density, plot.type=plot.type, ...)
  if (plot.type=='line') {
    abline(h=x$level, lty=2)
    abline(v=c(x$zeros), lty=2)
    for (i in 1:nrow(x$zeros)) {
      zero1 <- x$zeros[i,1]
      zero2 <- x$zeros[i,2]
      inside <- x$density$x >= zero1 & x$density$x <= zero2
      polygon(x=c(zero2, zero1, x$density$x[inside], zero2), y=c(0,0,x$density$y[inside], 0), density = polygon.control$density, angle = polygon.control$angle, border = polygon.control$border, col = polygon.control$col, lty = polygon.control$lty, fillOddEven = polygon.control$fillOddEven)
    }
  } else {
    warning('Not Yet Implemented for plot.type=circle')
  }
}

#############################################################
#
#	Internal functions for modal.region.circular
#       GNU General Public Licence 2.0 
#	Author: Claudio Agostinelli
#	E-mail: claudio@unive.it
#	Date: July, 5, 2011
#	Version: 0.5-1
#
#	Copyright (C) 2011 Claudio Agostinelli
#
#############################################################

### product of subsequent observations (like function diff.default)
prodseq <- function (x, lag = 1L, differences = 1L, ...) {
  ismat <- is.matrix(x)
  xlen <- if (ismat) 
    dim(x)[1L]
  else length(x)
    if (length(lag) > 1L || length(differences) > 1L || lag < 1L || differences < 1L) 
      stop("'lag' and 'differences' must be integers >= 1")
  if (lag * differences >= xlen) 
    return(x[0])
  r <- unclass(x)
  i1 <- -seq_len(lag)
  if (ismat) 
    for (i in seq_len(differences)) r <- r[i1, , drop = FALSE] * r[-nrow(r):-(nrow(r) - lag + 1L), , drop = FALSE]
  else for (i in seq_len(differences))
    r <- r[i1] * r[-length(r):-(length(r) - lag + 1L)]
  class(r) <- oldClass(x)
  return(r)
}

### Search for all crosses
allcrosses <- function(l, object, grid=seq(0,2*pi+0.01,0.01), tol=.Machine$double.eps^0.25, ...) {
  den <- approxfun(x=object$x, y=object$y)
  int <- function(x, l) {
    den(x) - l
  }
  x <- object$x
  y <- int(x=grid, l=l)
  sy <- sign(y)
  psy <- prodseq(sy)
  pos <- which(psy==-1)
  if (length(pos)) {
    intervals <- matrix(c(grid[pos], grid[pos+1]), ncol=2, byrow=FALSE)
    sy <- matrix(c(sy[pos], sy[pos+1]), ncol=2, byrow=FALSE)
    zeros <- rep(NA, nrow(intervals))
    for (i in 1:nrow(intervals)) {
      zeros[i] <- cross(l=l, object=object, lower=intervals[i,1], upper=intervals[i,2], tol=tol, ...)
    }
   if (length(zeros)%%2) {
      if (isTRUE(all.equal(zeros[1L]%%(2*pi), zeros[length(zeros)]%%(2*pi), tol=tol^0.9, scale=1)))
        zeros <- zeros[-length(zeros)]    
    }
    if (!length(zeros)%%2) {
      tsy <- apply(sy,2,prodseq)
      if (all(tsy==-1)) {
        if (sy[1,1]==1) {
          zeros <- c(0, zeros, 2*pi)
        }
        zeros <- matrix(zeros, ncol=2, byrow=TRUE)        
      } else {
        warning('Probably one zeros is missed, the zeros found are not in order')
      }
    } else {
      warning('The number of zeros is odd. At least one zero is missed')
    }
  } else {
    zeros <- NA            
  }
  return(zeros)
}

## Search the precise position of one cross
cross <- function(l, object, lower, upper, ...) {
#l: level
#object: an object from density.circular with results in radians
  den <- approxfun(x=object$x, y=object$y)
  int <- function(x, l) {
    den(x) - l
  }
  zero <- uniroot(int, l=l, lower=lower, upper=upper, ...)$root
  return(zero)
}

## Calculate areas under several disjoint intervals
area <- function(x, object, ...) {
#x: is a matrix with two columns
#object: an object from density.circular
#...: values passed to integrate function
  den <- approxfun(x=c(object$x-2*pi,object$x,object$x+2*pi), y=rep(object$y, 3))
  int <- function(x) integrate(f=den, lower=x[1], upper=x[2], ...)$value
  areas <- apply(x, 1, int)
  tot <- sum(areas)
  result <- list(tot=tot, areas=areas)
  return(result)
}


