
##==============================================================================
## Attaches a property to a 2D grid
##==============================================================================

setup.prop.2D <- function(func = NULL, value = NULL,  grid, 
                          y.func = func, y.value = value,...) {

  ## check input
  gn <- names(grid)
  if (! "x.mid"  %in% gn)
    stop("error in setup.2D.prop: grid should be a list that contains x.mid")
  if (! "x.int"  %in% gn)
    stop("error in setup.2D.prop: grid should be a list that contains x.int")
  if (! "y.mid"  %in% gn)
    stop("error in setup.2D.prop: grid should be a list that contains x.mid")
  if (! "y.int"  %in% gn)
    stop("error in setup.2D.prop: grid should be a list that contains x.int")

  if (is.null(func) && is.null(value))
    stop("error in setup.prop.2D: 'func' and 'value' should not both be NULL")

  if (is.null(y.func) && is.null(y.value))
    stop("error in setup.prop.2D: 'y.func' and 'y.value' should not both be NULL")

  if (!is.null(func) && !is.null(value))
    stop("error in setup.prop.2D: either 'func' or 'value' should be specified, not both")

  if (!is.null(y.func) && !is.null(y.value))
    stop("error in setup.prop.2D: either 'y.func' or 'y.value' should be specified, not both")

  Nx <- length(grid$x.mid)
  Ny <- length(grid$y.mid)  

  if (!is.null(value)) { # profile specification via constant value
    x.int <- matrix(nrow=Nx+1,ncol=Ny  ,data=value)
    y.int <- matrix(nrow=Nx  ,ncol=Ny+1,data=y.value)
    x.mid <- matrix(nrow=Nx  ,ncol=Ny  ,data=value)
    y.mid <- matrix(nrow=Nx  ,ncol=Ny  ,data=y.value)
  }

  if (!is.null(func)) { # profile specification via function
    if (is.vector(grid$x.mid)) {
    x.int <- outer(X=grid$x.int,Y=grid$y.mid, FUN=func, ...)
    y.int <- outer(X=grid$x.mid,Y=grid$y.int, FUN=y.func, ...)
    x.mid   <- outer(X=grid$x.mid,Y=grid$y.mid, FUN=func, ...)
    y.mid   <- outer(X=grid$x.mid,Y=grid$y.mid, FUN=y.func, ...)
    } else 
      stop("not yet implemented for matrix grid$x.mid")
   }
  
  Res <- list(x.mid = x.mid,
              y.mid = y.mid,  
              x.int = x.int,
              y.int = y.int)
  class(Res) <- "prop.2D"
  return(Res)
}

##==============================================================================
## S3 method: Plotting of a two-dimensional grid property
##==============================================================================

contour.prop.2D <- function(x, grid, xyswap = FALSE, filled = FALSE, ...) {
  if (! filled) {
  if (xyswap)
    contour(x=grid$y.mid,y=rev(-grid$x.mid),z=t(x$x.mid),...)
  else
    contour(x=grid$x.mid,y=grid$y.mid,z=x$x.mid,...)
  } else {
    if (xyswap)
      filled.contour(x=grid$y.mid,y=rev(-grid$x.mid),z=t(x$x.mid),...)
    else
      filled.contour(x=grid$x.mid,y=grid$y.mid,z=x$x.mid,...)
  }
}
