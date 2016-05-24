#' @title Plot arrow between variables
#' 
#' @description
#' Use this function to draw connecting arrows between manifest 
#' and latent variables.
#'
#' @param from An object of either class \code{"manifest"} or \code{"latent"}. 
#' This object is the origin of the arrow. 
#' @param to An object of either class \code{"manifest"} or \code{"latent"}. 
#' This object is the destination of the arrow.
#' @param start string to specify the starting direction of the arrow. 
#' Options are \code{"north"}, \code{"south"}, \code{"east"}, \code{"west"}.
#' @param end string to specify the ending direction of the arrow. Options
#' are \code{"north"}, \code{"south"}, \code{"east"}, \code{"west"}.
#' @param length of the edges of the arrow head (in inches).
#' @param angle from the shaft of the arrow to the edge of the arrow head.
#' @param code integer code, determining the kind of arrows to be drawn.
#' @param col color of the arrow.
#' @param lwd width of the arrow.
#' @param ... other arguments passed on to \code{\link{arrows}}.
#' @author Gaston Sanchez
#' @seealso \code{\link{manifest}}, \code{\link{latent}}, \code{\link{draw}}
#' @export
#' @examples
#'
#'  \dontrun{
#'  # latent variables
#'  attack = latent("ATTACK", x=0.35, y=0.7, rx=0.08, ry=0.06)
#'  defense = latent("DEFENSE", x=0.35, y=0.3, rx=0.08, ry=0.06)
#'  success = latent("SUCCESS", x=0.65, y=0.5, rx=0.08, ry=0.06)
#'
#'  # open wall
#'  wall()
#'  # draw latent variables
#'  draw(attack)
#'  draw(defense)
#'  draw(success)
#'  # add arrows
#'  arrow(from=attack, to=success, start="east", end="west")
#'  arrow(from=defense, to=success, start="east", end="west")
#'  }
#'
arrow <- 
function(from, to, start = "east", end = "west", length = 0.1,
         angle = 10, code = 2, col = "#d2def1", lwd = 3, ...)
{
  ini = grep(start, c("north", "south", "east", "west")) + 5
  fin = grep(end, c("north", "south", "east", "west")) + 5
  x0 = from[[ini]][1]
  y0 = from[[ini]][2]  
  x1 = to[[fin]][1]
  y1 = to[[fin]][2]
  arrows(x0=x0, y0=y0, x1=x1, y1=y1, length=length, angle=angle, 
         code=code, col=col, lwd=lwd, ...)
}
