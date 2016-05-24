#' Create a SpatialPointsDataFrame describing Metapop population centroids
#' 
#' Create a \code{SpatialPointsDataFrame} representing the centroid of each
#' population, with attributes: \code{pop} (population name), \code{time} (the
#' time step), and \code{N} (the mean population size).
#' 
#' @param mp A RAMAS Metapop .mp file containing simulation results.
#' @param coords An object containing population coordinates. This object can be
#'   created by using \code{\link{mp2xy}}
#' @param start The value of the first timestep. If timesteps are not in 
#'   increments of 1, it may be best to use \code{start=1}, in which case 'time'
#'   in the resulting shapefile's attribute table will refer to the timestep 
#'   number.
#' @param s_p4s (Optional) The coordinate reference system of the source
#'   cordinates given in \code{coords}. These can be supplied as a \code{CRS} 
#'   object or as a proj4 string.
#' @param t_p4s (Optional) The target coordinate reference system to which
#'   coordinates will be projected, if supplied. These can be supplied as a
#'   \code{CRS} object or as a proj4 string.
#' @return A \code{SpatialPointsDataFrame}.
#' @keywords spatial
#' @seealso \code{\link{mp2xy}}
#' @importFrom sp coordinates CRS proj4string spTransform
#' @importFrom raster crs
#' @export
#' @examples
#' mp <- system.file('example.mp', package='mptools')
#' r <- system.file('example_001.tif', package='mptools')
#' coords <- mp2xy(mp, r, 9.975)
#' spdf <- mp2sp(mp, coords, start=2000)
#' spdf
mp2sp <- function(mp, coords, start, s_p4s, t_p4s) {
  res <- results(mp)
  sites <- coords[, c('pop', 'x', 'y')]
  N <- as.data.frame(t(res$results[, 'mean', -1]))
  names(N) <- start + seq_along(N) - 1
  if(!identical(row.names(N), sites$pop)) 
    stop('Something went wrong. Please contact the package maintainer.')
  N <- utils::stack(N)[, 2:1]
  names(N) <- c('time', 'N')
  spdf <- cbind(sites, N)
  sp::coordinates(spdf) <- ~x+y
  if(!missing('s_p4s')) {
    tryCatch(raster::crs(s_p4s), error=function(e) 
      stop('proj4string not recognised: ', s_p4s, call.=FALSE))
    sp::proj4string(spdf) <- s_p4s 
  }
  if(!missing('t_p4s')) {
    if(missing('s_p4s'))
       stop('If t_p4s is supplied, s_p4s must also be supplied', call.=FALSE)
    tryCatch(raster::crs(t_p4s), error=function(e) 
      stop('proj4string not recognised: ', t_p4s, call.=FALSE))
    spdf <- sp::spTransform(spdf, sp::CRS(t_p4s)) 
  }
  spdf
}
