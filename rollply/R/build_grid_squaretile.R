#'
#' @title Create a grid with square tiles
#'
#' @description Create a rectangular grid of regularly-spaced points (square
#'             tiles).
#' 
#' @param coords A matrix or \code{data.frame} of coordinates with two columns
#' @param npts The approximate total number of points of the output grid
#' @param pad Padding on each dimension (a positive number makes a grid
#'            that is larger than the ranges of the coordinates).
#' @param ... other arguments are silently ignored
#'
#' @return The coordinates of a grid of points as a \code{data.frame} with
#'         approximately \code{npts} rows and \code{ncol(coords)} columns. Names
#'         are transfered from the \code{coords} data frame.
#' 
#' @details This function creates a grid that covers a set of points. The 
#'          distance between points is the same on all dimensions (tiles are 
#'          squared. It is only implemented for 2D grid so \code{coords} must
#'          have at most two columns.
#' 
#' @family grid building functions
#'
#' @export
build_grid_squaretile <- function(coords, npts, pad = 0, ...) {

  build_grid_check_vars(coords, npts)

  coords.ranges <- apply(coords, 2, range)
  ndims <- ncol(coords)

  if (ncol(coords)!=2) {
    stop('This type of grid is only implemented for 2 dimensions.')
  }

  # Get the number of points to add on each dimension
  npts.dims <- get_npts_dims(coords.ranges, npts, pad)

  # Get the size of one tile
  tile.sizes <- apply(coords.ranges, 2, diff) / npts.dims
  tile.size <- mean(tile.sizes) # this is our tile size


  # Build grid seed to feed to expand.grid
  grid.seed <- list()
  for (col in seq.int(ncol(coords))) {
    grid.seed[[col]] <- seq(min(coords.ranges[ ,col]) - pad,
                            max(coords.ranges[ ,col]) + pad,
                            by = tile.size)
  }

  names(grid.seed) <- colnames(coords)

  expand.grid(grid.seed, KEEP.OUT.ATTRS = FALSE)
}

get_npts_dims <- function(coords.ranges, npts, pad) {
  # Given npts to put in a L*l rectangle of size ratio r (=L/l), then
  # m=sqrt(r*npts) and n=r*sqrt(npts) to have m x n tiles.

  sides.l <- apply(coords.ranges, 2, diff) + pad
  r <- max(sides.l) / min(sides.l)
  npts.dims <- c(sqrt(npts/r), sqrt(npts*r))[order(sides.l)]
  npts.dims <- round(npts.dims)

  return(npts.dims)
}

