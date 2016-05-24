
##==============================================================================
# Creation of a two-dimensional finite difference grid
##==============================================================================

setup.grid.2D <- function(x.grid=NULL, y.grid=NULL) {

## check input
  gn <- names(x.grid)
  if (! "x.up"   %in% gn)
    stop("error in setup.2Dgrid: x.grid should be a list that contains x.up")
  if (! "x.down" %in% gn)
    stop("error in setup.2Dgrid: x.grid  should be a list that contains x.down")
  if (! "x.mid"  %in% gn)
    stop("error in setup.2Dgrid: x.grid  should be a list that contains x.mid")
  if (! "x.int"  %in% gn)
    stop("error in setup.2Dgrid: x.grid  should be a list that contains x.int")
  if (! "dx"     %in% gn)
    stop("error in setup.2Dgrid: x.grid  should be a list that contains dx")
  if (! "dx.aux" %in% gn)
    stop("error in setup.2Dgrid: x.grid  should be a list that contains dx.aux")
  if (! "N"      %in% gn)
    stop("error in setup.2Dgrid: x.grid  should be a list that contains N")

  gn <- names(y.grid)
  if (! "x.up"   %in% gn)
    stop("error in setup.2Dgrid: y.grid should be a list that contains x.up")
  if (! "x.down" %in% gn)
    stop("error in setup.2Dgrid: y.grid  should be a list that contains x.down")
  if (! "x.mid"  %in% gn)
    stop("error in setup.2Dgrid: y.grid  should be a list that contains x.mid")
  if (! "x.int"  %in% gn)
    stop("error in setup.2Dgrid: y.grid  should be a list that contains x.int")
  if (! "dx"     %in% gn)
    stop("error in setup.2Dgrid: y.grid  should be a list that contains dx")
  if (! "dx.aux" %in% gn)
    stop("error in setup.2Dgrid: y.grid  should be a list that contains dx.aux")
  if (! "N"      %in% gn)
    stop("error in setup.2Dgrid: y.grid  should be a list that contains N")

## Packaging of results
  Nx <- length(x.grid$x.mid)
  Ny <- length(y.grid$x.mid)
  
  Res <- list(x.up = x.grid$x.up,
				    x.down = x.grid$x.down,
				    x.mid = x.grid$x.mid,                               # position of centre of the grid cells, vector of length N
            x.int  = x.grid$x.int,                              # position grid cell interfaces , vector of length N+1
            dx = matrix(nrow=Nx, ncol=Ny, x.grid$dx),           # thickness of grid cells , vector length N
            dx.aux = matrix(nrow=Nx+1, ncol=Ny, x.grid$dx.aux), # auxiliary vector with distances between centre of adjacent cells, first and last: half of cell size, vector of length N+1
            x.N = x.grid$N,         # number of vertical grid layers
            y.up = y.grid$x.up,
				    y.down = y.grid$x.down,
				    y.mid = y.grid$x.mid,   # position of centre of the grid cells, vector of length N
            y.int  = y.grid$x.int,  # position grid cell interfaces , vector of length N+1
            dy = matrix(nrow=Nx, ncol=Ny, data=y.grid$dx, byrow=TRUE),        # thickness of grid cells , vector length N
            dy.aux = matrix(nrow=Nx, ncol=Ny+1, y.grid$dx.aux, byrow=TRUE), # auxiliary vector with distances between centre of adjacent cells, first and last: half of cell size, vector of length N+1
            y.N = y.grid$N)         # number of horizontal grid layers

  class(Res) <- "grid.2D"
  return(Res)
}
