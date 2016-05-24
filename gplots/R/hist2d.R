# $Id: hist2d.R 1555 2012-06-07 16:08:18Z warnes $

hist2d <- function(x,
                   y=NULL,
                   nbins=200,
                   same.scale=FALSE,
                   na.rm=TRUE,
                   show=TRUE,
                   col=c("black", heat.colors(12)),
                   FUN=base::length,
                   xlab,
                   ylab,
                   ... )
  {
    if(is.null(y))
      {
        if(ncol(x) != 2) stop("If y is ommitted, x must be a 2 column matirx")
        y <- x[,2]
        x <- x[,1]
      }

    if(length(nbins)==1)
      nbins <- rep(nbins,2)

    nas <- is.na(x) | is.na(y)

    if(na.rm)
      {
        x <- x[!nas]
        y <- y[!nas]
      }
    else
      stop("missinig values not permitted if na.rm=FALSE")

    if(same.scale)
      {
        x.cuts <- seq( from=min(x,y), to=max(x,y), length=nbins[1]+1)
        y.cuts <- seq( from=min(x,y), to=max(x,y), length=nbins[2]+1)
      }
    else
      {
        x.cuts <- seq( from=min(x), to=max(x), length=nbins[1]+1)
        y.cuts <- seq( from=min(y), to=max(y), length=nbins[2]+1)
      }

    index.x <- cut( x, x.cuts, include.lowest=TRUE)
    index.y <- cut( y, y.cuts, include.lowest=TRUE)

    ## tapply is faster than old for() loop, and allows
    ## use of any user-specified summary function 
    m <- tapply(x,list(index.x,index.y),FUN)

    ## If we're using length, set empty cells to 0 instead of NA
    if(identical(FUN,base::length))
      m[is.na(m)] <- 0

    if(missing(xlab)) xlab <- deparse(substitute(xlab))
    if(missing(ylab)) ylab <- deparse(substitute(ylab))
    
    if(show)
      image( x.cuts, y.cuts, m, col=col, xlab=xlab, ylab=ylab, ...)

    midpoints <- function(x) (x[-1]+x[-length(x)])/2
    
    retval <- list()
    retval$counts <- m
    retval$x.breaks = x.cuts
    retval$y.breaks = y.cuts
    retval$x = midpoints(x.cuts)
    retval$y = midpoints(y.cuts)
    retval$nobs = length(x)
    retval$call <- match.call()
    class(retval) <- "hist2d"
    retval
  }
