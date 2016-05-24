#'
#' @title Regular grid in alpha-hull
#'
#' @description Create a grid within the alpha hull of a set of points.
#'
#' @details This function computes the alpha hull of a set of points, then
#'   iteratively finds the best grid of \code{npts} points fitting in the hull.
#' 
#' @param coords A matrix or data.frame of coordinates with two columns
#' 
#' @param npts The approximate number of points of the requested grid
#' 
#' @param pad ignored argument
#' 
#' @param grid_opts A list with component \code{alpha}, \code{error_tol}, 
#'   \code{run_max} and \code{verbose} that controls the
#'   computation of the grid.
#' 
#' @return The coordinates of a grid of points as a \code{data.frame} with
#'   approximately \code{npts} rows and \code{ncol(coords)} columns. Names
#'   are transfered from the \code{coords} data frame.
#' 
#' @details This function computes the alpha hull of a set of points which is 
#'   an estimate of its overall outline. The intricateness of this outline 
#'   is determined by the \code{alpha} parameter: a value close to zero will 
#'   produce a very intricate outline and a high value will produce a coarse 
#'   outline (similar to the convex hull). It is highly recommended to try 
#'   different values of \code{alpha} to see what is most appropriate for one's
#'   datasets. 
#'   
#'   Parameter \code{error_tol} controls the acceptable relative error in the 
#'   grid number of points (e.g. 0.1 for 10%). \code{run_max} and \code{verbose}
#'   control the iteration process.
#'   
#' @seealso \code{\link[alphahull]{ahull}}
#' 
#' @family grid building functions
#' 
#' @examples 
#' 
#' data(meadow)
#' 
#' # Plot a grid with diagnostics
#' grid <- build_grid_ahull_fill(meadow[ ,c('lon','lat')], 10000,
#'                               grid_opts = list(verbose = TRUE))
#' plot(grid, pch = 21)
#' 
#' # See the influence of alpha parameter
#' par(mfrow = c(1, 4))
#' for ( alpha in c( .3, 1, 4, 10) ) {
#'   grid <- build_grid_ahull_fill(meadow[ ,c('lon','lat')], 10000,
#'                                 grid_opts = list(alpha = alpha))
#'   plot(grid, pch = 21)
#'   title(paste0("alpha=",alpha))
#' }
#' 
#' @export
build_grid_ahull_fill <- function(coords, npts,
                                  pad = NULL, # pad is ignored for this function
                                  grid_opts = list(alpha = .3,
                                                   error_tol = .05,
                                                   run_max = 20,
                                                   verbose = FALSE)) {

  if ( !requireNamespace("alphahull", quietly = TRUE) ) {
    stop('Alpha-hull-based grids require package alphahull')
  }

  build_grid_check_vars(coords, npts)

  if (ncol(coords)!=2) {
    stop('Incorrect number of columns in the coordinates provided (required 2 -',
         paste(" got:", ncol(coords)), ")")
  }

  # Remove duplicates otherwise ahull throws an error
  coords <- coords[ ! duplicated(coords), ]
  
  # Rescale if dimensions have very big numbers. This fixes a bug in 
  #   inahull at the cost of being less precise
  scale_factor <- max(apply(coords, 2, function(X) diff(range(X))))
  coords <- coords / scale_factor
  
  # Take parameters into account
  opts <- list(alpha = .3, error_tol = .05, run_max = 20, verbose = FALSE)
  opts[names(grid_opts)] <- grid_opts # alter defaults
  
  
  # We build an alpha hull of our x/y points.
  if ( opts[['verbose']] ) cat('Building alphahull...\n')
  coords.hull <- alphahull::ahull(coords, alpha = opts[["alpha"]]) # opts[['alpha']])
  
  # Iterate to find the closest number of points
  error <- -1
  run <- 1
  npts_target <- npts
  to_keep <- TRUE
  while ( error < - opts[['error_tol']] && 
          run <= opts[['run_max']] && 
          any(to_keep) ) {
    
    # Compute a grid with the given number of points
    rect_grid <- build_grid_squaretile(coords, npts)
    
    # Strip points and see how many fall in alphahull
    to_keep <- inahull_cpp_multiple(coords.hull, 
                                    rect_grid[[1]], 
                                    rect_grid[[2]])
    
    error <- (sum(to_keep) - npts_target) / npts_target
    
    if ( opts[['verbose']] ) {
      cat(paste0('run ', run,', error=', round(error, digits = 2),
                     '% (', sum(to_keep), ' points)\n'));
    }
    
    # Increase number of points
    npts <- round(npts - error * npts) # error is < 0
    run <- run + 1
  }
  
  if ( !any(to_keep) ) { 
    warning("No points of the grid fell in the alphahull: maybe parameter 
             alpha is too small ?")
  }
  
  # Scale back to original size and return grid
  grid <- rect_grid[to_keep, ] * scale_factor
  
  return(grid)
}
