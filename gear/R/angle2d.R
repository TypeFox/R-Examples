#' Determine angle
#' 
#' \code{angle2d} determines the angle between pairs of coordinates in degrees or radians.  The coordinats are assumed to be in 2d space.
#' 
#' Note that the angle is between the actual pairs of points, not the angle between the vectors extending from the origin to the points.  e.g., the angle between cbind(0, 1) and cbind(1, 1) would be 90 degrees, not 45.  Sign of the direction not accounted for, e.g., a -135 degree angle is rotated by 180 degrees to become a 45 degree angle.  All angles returned are in the interval [0, 180].
#' 
#' @param coords1 An \eqn{N \times 2} matrix of spatial coordinates.
#' @param coords2 An \eqn{N \times 2} matrix of spatial coordinates.
#' @param radians A logical indicating whether degrees or radians should be returned.  Default is FALSE, meaning return angle in degrees.
#' @return Returns a vector of angles. 
#' @author Joshua French
#' @export
#' @examples 
#' coords1 = matrix(0, nrow = 8, ncol = 2)
#' coords2 = cbind(c(2, 2, 0, -2, -2, -2, 0, 2), c(0, 2, 2, 2, 0, -2, -2, -2))
#' angle2d(coords1, coords2)
#' angle2d(coords1, coords2, radians = TRUE)

angle2d <- function(coords1, coords2, radians = FALSE){
  # sanity checking of arguments
  if(ncol(coords1) != 2) stop(ncol(coords1) != 2)
  if(ncol(coords2) != 2) stop(ncol(coords2) != 2)
  if(nrow(coords1) != nrow(coords2)) 
    stop("nrow(coords1) != nrow(coords2)")
  if(!is.logical(radians) || length(radians) != 1) 
    stop("radians must be a single logical value")
  # determine directional y- and x-vectors (moving from 1st coordinate to 2nd)
  y = coords2[,2] - coords1[,2] 
  x = coords2[,1] - coords1[,1]
  # determine the angle between the coordinates, based on these vectors
  ang = atan(y/x) 
  # since we don't care about the sign of the angle (e.g., -45 degrees is equal to 135 degrees)
  # convert negative radians to positive by rotating 180 degrees
  ang[ang < 0] = ang[ang < 0] + pi
  # convert to degrees if necessary
  ang * ifelse(!radians, 180/pi, 1)
}
