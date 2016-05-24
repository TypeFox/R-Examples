
##==============================================================================
## Attaches a property to a 1D grid
##==============================================================================

setup.prop.1D <- function(func = NULL, value = NULL, xy = NULL,
  interpolate = "spline", grid, ...) {

## check input
  gn <- names(grid)
  if (! "x.up"   %in% gn)
    stop("error in setup.prop.1D: grid should be a list that contains x.up")
  if (! "x.down" %in% gn)
    stop("error in setup.prop.1D: grid should be a list that contains x.down")
  if (! "x.mid"  %in% gn)
    stop("error in setup.prop.1D: grid should be a list that contains x.mid")
  if (! "x.int"  %in% gn)
    stop("error in setup.prop.1D: grid should be a list that contains x.int")
  if (! "dx"     %in% gn)
    stop("error in setup.prop.1D: grid should be a list that contains dx")
  if (! "dx.aux" %in% gn)
    stop("error in setup.prop.1D: grid should be a list that contains dx.aux")
  if (! "N"      %in% gn)
    stop("error in setup.prop.1D: grid should be a list that contains N")


  if (is.null(xy) && is.null(value) && is.null(func))
    stop("error in setup.prop.1D: function, value and xy cannot be NULL together")

## profile specification via function
  if (! is.null(func)) {   
  	y.int <- func(grid$x.int,...)
    ## KS ##
  	if (length(y.int)==0)
      stop (" error in setup.prop.1D: 'func' should return a vector ")
      
	  y.mid <- func(grid$x.mid,...)
  }

  if (! is.null(value)) { # profile specification via value
	  y.int <- rep(value,times=length(grid$x.int))
	  y.mid <- rep(value,times=length(grid$x.mid))
  }

## profile speficication via data series input
  if (! is.null(xy)) { 
    if (! is.matrix(xy))
      stop("error in setup.prop.1D: xy should be a 2-columned matrix or NULL")
	  if (!(interpolate %in% c("linear","spline")))
      stop("error in setup.prop.1D: 'interpolate' not properly specified, should be 'linear' or 'spline'")
		
    if (interpolate=="linear") {
      # check the range of input data-series (inr))
      # this range should span the range of grid$x.int (outr)
      outr <- range(grid$x.int) # output range
      inr  <- range(xy[,1])     # input range
	    # extend the xy matrix if necessary
	    if (outr[1]<inr[1]) {
        xy <- rbind(c(outr[1],xy[which.min(xy[,1]),2]),xy)
      }
      if (outr[2]>inr[2]) {
        xy <- rbind(xy,c(outr[2],xy[which.max(xy[,1]),2]))
      }
      # linear interpolation
    	y.int <- approx(x=xy[,1],y=xy[,2],xout = grid$x.int)$y
      y.mid <- approx(x=xy[,1],y=xy[,2],xout = grid$x.mid)$y
    } else {

	    # spline interpolation
	    f <- splinefun(xy[,1],xy[,2])
      y.int <- f(grid$x.int)
      y.mid <- f(grid$x.mid)
    }
  }

  Res <- list(mid  = y.mid,
              int  = y.int)
  class(Res) <- "prop.1D"
  return(Res)
}


##==============================================================================
## S3 method: Plotting of a one-dimensional grid property
##==============================================================================
setdots <- function (dots, default) {
  dots <- if (is.null(dots)) default else dots
}


plot.prop.1D <- function(x, grid, xyswap =FALSE, ...) {
  dots <- list(...)
     
  if (xyswap) {
    dots$ylim <- setdots(dots$ylim, rev(range(grid$x.int)))
    dots$xlab <- setdots(dots$xlab, "prop")
    dots$ylab <- setdots(dots$ylab, "x")
    List <- alist(x = x$int, y = grid$x.int)
  } else {
    dots$xlab <- setdots(dots$xlab, "x")
    dots$ylab <- setdots(dots$ylab, "prop")
    List <- alist(x = grid$x.int, y = x$int)
  }    
  do.call(plot, c(List, dots))

}
