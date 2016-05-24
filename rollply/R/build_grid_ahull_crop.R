#' @title Regular grid in alpha-hull
#' 
#' @description Create a grid within the alpha hull of a set of points.
#' 
#' @details This function creates a rectangular grid over a set of points, then 
#'   computes the alpha hull of the set of points, and discards all the points 
#'   of the new grid that fall outside of the hull.
#' 
#' @param coords A matrix or data.frame of coordinates with two columns
#' 
#' @param npts The number of points before cropping to the alpha hull shape
#' 
#' @param grid_opts A list with component \code{alpha} that controls the shape 
#'   of the alpha hull (see details).
#'  
#' @param pad Ignored
#' 
#' @return The coordinates of a grid of points as a \code{data.frame} with
#'         \code{ncol(coords)} columns. Names are transfered from the
#'         \code{coords} data frame.
#' 
#' @details This function computes the alpha hull of a set of points which is 
#'   an estimate of its overall outline. The intricateness of this outline 
#'   is determined by the \code{alpha} parameter: a value close to zero will 
#'   produce a very intricate outline and a high value will produce a coarse 
#'   outline (similar to the convex hull). It is highly recommended to try 
#'   different values of \code{alpha} to see what is most appropriate for one's
#'   datasets. 
#' 
#' @seealso \code{\link[alphahull]{ahull}}
#' 
#' @family grid building functions
#'
#'@export
build_grid_ahull_crop <- function(coords, npts,
                                  pad = NULL, # pad is ignored for this function
                                  grid_opts = list(alpha = .3, 
                                                   scale_alpha = TRUE)) {

  build_grid_ahull_fill(coords, npts,
                        grid_opts = c(grid_opts,
                                      list(run_max = 1,
                                           verbose = FALSE)))
}
