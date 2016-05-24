#' Distance between points
#' 
#' Calculate distance between points on planar surface
#' 
#' @details The function is quite simple: \code{sqrt((xref - xpt)^2 + (yref - ypt)^2)}
#' 
#' @return vector with the distances
#' @author Berry Boessenkool, \email{berry-b@@gmx.de}, 2012
#' @seealso \code{\link[spatstat]{nndist}} in the package \code{spatstat} for distance to nearest neighbour
#' @keywords spatial
#' @export
#' @examples
#' 
#' A <- c(3,  9,-1)
#' B <- c(7, -2, 4)
#' plot(A,B)
#' text(A,B, paste0("P",1:3), adj=1.1)
#' points(3,5, col=2, pch=16)
#' segments(3,5, A,B)
#' distance(A,B, 3,5)
#' text(c(3.2,6,1), c(6,1,4), round(distance(A,B, 3,5),2) )
#' 
#' @param xpt vector with x-coordinate(s) of point(s)
#' @param ypt ditto for y
#' @param xref single x coordinate of reference point
#' @param yref ditto for y
#' 
distance <- function(
   xpt,
   ypt,
   xref,
   yref)
{
sqrt((xref-xpt)^2 + (yref-ypt)^2)
}
