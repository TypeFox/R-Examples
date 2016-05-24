#' Determine whether circles intersection
#' 
#' \code{circles.intersect} determines whether circles intersect with each other.
#' 
#' The algorithm is based on the premise that two circles intersect if, and only if, the distance between their centroids is between the sum and the difference of their radii.  I have squared the respective parts of the inequality in the implemented algorithm.
#' 
#' @param coords A matrix of coordinates containined the centroid of each circle.
#' @param r A vector containing the radii of the circles.  The have length equal to the number of rows of
#' \code{coords}.
#' @return Returns a matrix of logical values indicating whether the circles intersect.
#' @author Joshua French
#' @importFrom SpatialTools dist1
#' @export
#' @examples 
#' # first two intersect, then next two, last doesn't intersect any other
#' co = cbind(c(1, 2, 5, 6, 9), c(1, 2, 5, 6, 9))
#' r = c(1.25, 1.25, 1.25, 1.25, 1.25)
#' circles.intersect(co, r)
#' # nested circles
#' co = matrix(rep(0, 4), nrow = 2)
#' r = c(1, 1.5)
#' circles.intersect(co, r)

circles.intersect <- function(coords, r)
{
  d = SpatialTools::dist1(coords)
  if(length(r) != nrow(coords)) stop("length(r) must be equal to nrow(coords)")
  rmat1 = matrix(r, nrow = length(r), ncol = length(r))
  rmat2 = t(rmat1)
  lb = (rmat1 - rmat2)^2
  ub = (rmat1 + rmat2)^2
  return(lb <= d^2 & d^2 <= ub)
}

