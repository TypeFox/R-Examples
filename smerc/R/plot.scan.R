#' Plots object of class \code{scan}. 
#' 
#' Plots clusters (the centroids of the regions in each cluster) in different colors.  The most likely cluster is plotted with solid red circles by default.  Points not in a cluster are black open circles.  The other cluster points are plotted with different symbols and colors.  
#'
#' @param x An object of class scan to be plotted.
#' @param ... Additional graphical parameters passed to \code{plot} function.
#' @param ccol Fill color of the plotted points.  Default is NULL, indicating red for the most likely cluster, and col = 3, 4, ..., up to the remaining number of clusters.
#' @param cpch Plotting character to use for points in each cluster.  Default is NULL, indicating pch = 20 for the most likely cluster and then pch = 2, 3, .., up to the remaining number of clusters.
#' @param usemap Logical indicating whether the maps::map function should be used to create a plot background for the coordinates.  Default is FALSE.  Use TRUE if you have longitude/latitude coordinates.
#' @param mapargs A list of arguments for the map function.
#' @importFrom graphics plot points
#' @import maps
#' @method plot scan
#' @export
#' @seealso \code{\link[maps]{map}}
#' @examples
#' data(nydf)
#' coords = with(nydf, cbind(longitude, latitude))
#' out = scan.test(coords = coords, cases = floor(nydf$cases), 
#'                 pop = nydf$pop, nsim = 49, 
#'                 lonlat = TRUE, alpha = 0.12,
#'                 parallel = FALSE) 
#' ## plot output for new york state
#' # specify desired argument values
#' mapargs = list(database = "state", region = "new york", 
#' xlim = range(out$coords[,1]), ylim = range(out$coords[,2]))
#' # needed for "state" database (unless you execute library(maps))
#' data(stateMapEnv, package = "maps") 
#' plot(out, usemap = TRUE, mapargs = mapargs)

plot.scan = function(x, ..., ccol = NULL, cpch = NULL, usemap = FALSE, mapargs = list())
{
  if(class(x) != "scan") stop("x should be an scan object from an appropriate function, e.g., scan.test")
  
  # number of centroids
  nc = length(x$clusters)
  
  # set default values
  if(is.null(ccol)) ccol = c(2:(nc + 1))
  if(is.null(cpch)) cpch = rep(20, nc)
  if(nc > 1) cpch[2:nc] = 3:(nc + 1)
  
  # more sanity checking
  if(length(ccol) != nc) stop("if specified, ccol must have length equal to length(x$clusters)")
  if(length(cpch) != nc) stop("if specified, cpch must have length equal to length(x$clusters)")
  
  # extract coordinates and cluster coordinates
  coords = x$coords
  # ccoords = matrix(0, nrow = nc, ncol = 2)
  # for(i in 1:nc) ccoords[i, ] = coords[x$clusters[[i]]$loc, ]
  
  if(usemap)
  { do.call(maps::map, mapargs) }else
  {
    plot(coords, ...)
  }
  
  graphics::points(coords, ...)
  # plot clusters
  for(i in 1:nc)
  {
    graphics::points(coords[x$clusters[[i]]$locids,1], coords[x$clusters[[i]]$locids,2], col = ccol[i], pch = cpch[i])
  }
}