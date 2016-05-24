#' Color clusters
#' 
#' \code{color.clusters} is a simple helper function that makes it easier to color clusters of regions produced by an appropriate method, e.g., \code{scan.test} or \code{uls.test}.  Regions/clusters that are not part of any cluster have no color.
#' 
#' @param x An object of class scan produced by a function such as \code{scan.test}.
#' @param col A vector of colors to color the clusters in \code{x}.  Should have same length as the number of clusters in \code{x}.
#' @return Returns a vector with colors for each region/centroid for the data set used to construct \code{x}.  
#' @author Joshua French
#' @export
#' @examples 
#' data(nydf)
#' coords = with(nydf, cbind(longitude, latitude))
#' out = scan.test(coords = coords, cases = floor(nydf$cases), 
#'                   pop = nydf$pop, alpha = 0.12, lonlat = TRUE,
#'                   nsim = 49)
#' data(nypoly)
#' library(sp)
#' plot(nypoly, col = color.clusters(out))

color.clusters = function(x, col = 2:(length(x$clusters) + 1))
{
  if(class(x) != "scan") stop("x should be an object of class scan.")
  if(length(x$clusters) != length(col)) stop("The number of colors must match the number of clusters.")
  
  mycol = numeric(nrow(x$coords))
  for(i in seq_along(x$clusters))
  {
    mycol[x$clusters[[i]]$loc] = col[i]
  }
  return(mycol)
}