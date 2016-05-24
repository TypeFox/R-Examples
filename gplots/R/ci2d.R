# $Id: ci2d.R 1471 2011-08-16 01:03:31Z warnes $

## first(...) selects the first element of which(...)
first <- function(x,...)
  {
  w <- which(x,...)
  if(length(x)>1)
    w[1]
  else
    w
}

## first(...) selects the first element of which(...)
last <- function(x,...)
  {
  w <- which(x,...)
  if(length(x)>1)
    rev(w)[1]
  else
    w
}

## non-parametric 2 dimensional approximate confidence interval
ci2d <- function(x,
                 y = NULL,
                 nbins=51,
                 method=c("bkde2D","hist2d"),
                 bandwidth,
                 factor=1.0,
                 
                 ci.levels=c(0.50,0.75,0.90,0.95,0.975),
                 
                 show=c("filled.contour","contour","image","none"),
                 col=topo.colors(length(breaks)-1),
                 show.points=FALSE,
                 pch=par("pch"),
                 points.col="red",
                 xlab, ylab, 
                 ...)
  {

    show <- match.arg(show)
    method <- match.arg(method)
    breaks <- unique(c(0, ci.levels, 1.0))

    # get labels for x and y
    if (missing(xlab)) 
      xlab <- if (missing(x)) 
        ""
      else deparse(substitute(x))
    if (missing(ylab)) 
      ylab <- if (missing(y)) 
        ""
      else deparse(substitute(y))

    if(!is.null(y))
      x <- cbind(x,y)

    if(method=="hist2d")
      {
        h2d <- hist2d(x,
                      show=FALSE,
                      nbins=nbins,
                      ...)
        ## normalize
        h2d$density <- h2d$counts / sum(h2d$counts, na.rm=TRUE)
      }
    else if (method=="bkde2D")
      {
        if(length(nbins)==1)
          nbins <- c(nbins, nbins)

        if(missing(bandwidth))
          {
            h.x = dpik(x[,1])
            h.y = dpik(x[,2])
            bandwidth <-  c(h.x, h.y)
          }

        est <- bkde2D(x,
                      bandwidth=bandwidth*factor,
                      gridsize=nbins,
                      ...
                      )

        h2d <- list()
        h2d$x <- est$x1
        h2d$y <- est$x2
        h2d$counts <- est$fhat
        h2d$nobs <- nrow(x)
        h2d$density <- est$fhat / sum(est$fhat) # normalize
      }
    else
      stop("Unknown method: '", method, "'")

    uniqueVals <- rev(unique(sort(h2d$density)))
    cumProbs <- sapply(uniqueVals,
                       function(val) sum( h2d$density[h2d$density>=val] ) )
    names(cumProbs) <- uniqueVals
    h2d$cumDensity <- matrix(nrow=nrow(h2d$density), ncol=ncol(h2d$density))
    h2d$cumDensity[] <- cumProbs[as.character(h2d$density)]

    if(show=="image")
      {
        image(h2d$x, h2d$y, h2d$cumDensity,
              xlab=xlab, ylab=ylab,
              breaks=breaks, col=col)
        if(show.points)
            points(x[,1], x[,2], pch=pch, col=points.col);          
      }
    else if(show=="filled.contour")
      {
        if(show.points)
          plot.title <- function() {
            points(x[,1], x[,2], pch=pch, col=points.col);
          }
        else
          plot.title <- function() {}

        
        filled.contour(h2d$x, h2d$y, h2d$cumDensity,
                       levels=breaks,
                       col=col,
                       xlab=xlab,
                       ylab=ylab,
                       plot.title=plot.title(),
                       key.title=title("\nCI Level"),
                       key.axes=axis(4, at=breaks)
                       )
       }
    else if(show=="contour")
      {
        tmpBreaks <- breaks[breaks<1]  # avoid having 1.0 line
        contour(h2d$x, h2d$y, h2d$cumDensity,
                levels=tmpBreaks,
                labels=tmpBreaks,
                xlab=xlab,
                ylab=ylab,
                nlevels=length(tmpBreaks),
                col=col
                )
        if(show.points)
            points(x[,1], x[,2], pch=pch, col=points.col);          
      }

    h2d$contours <- contourLines(h2d$x, h2d$y, h2d$cumDensity,
                                 levels=breaks, nlevels=length(breaks))

    # use the confidence level value as the name in the contour list
    names(h2d$contours) <- sapply(h2d$contours, function(x) x$level)

    # convert each contour into a (x,y) dataframe 
    h2d$contours <- lapply( h2d$contours,
                            function(J) data.frame(x=J$x, y=J$y) )

    h2d$call <- match.call()
    class(h2d) <- "ci2d"
    
    invisible(h2d)
  }

