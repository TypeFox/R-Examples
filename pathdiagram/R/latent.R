#' @title Set specifications of a latent variable
#' 
#' @description
#' Use this function to specify the graphic characteristics of a latent variable.
#' The specifications will be used by the function \code{draw} to plot latent
#' variables (in a path diagram).
#' 
#' @details
#' Latent variables are drawn as ellipses using the function 
#' \code{\link{plotellipse}}
#'
#' @param label A character string with the label to be displayed.
#' @param x x-axis coordinate for center of ellipse.
#' @param y y-axis coordinate for center of ellipse.
#' @param rx long radius of ellipse.
#' @param ry short radius of ellipse.
#' @param border color of the border.
#' @param lwd width of border line.
#' @param fill color to fill the ellipse
#' @param col color of the label.
#' @param cex numeric character expansion of the label. 
#' @param vfont font family of the label.
#' @param font An integer specifying which font to use for the label.
#' See \code{\link{par}}
#' @param family The name of a font family for drawing text. 
#' Standard values are \code{"serif"}, \code{"sans"} and \code{"mono"}.
#' @return An object of class \code{"latent"}, which is a list with the
#' specified parameters to draw latent variables.
#' @author Gaston Sanchez
#' @seealso \code{\link{manifest}}, \code{\link{draw}}
#' @export
#' @examples
#'
#'  \dontrun{
#'  # latent variables
#'  attack = latent("ATTACK", x=0.35, y=0.7, rx=0.08, ry=0.06)
#'  defense = latent("DEFENSE", x=0.35, y=0.3, rx=0.08, ry=0.06)
#'  success = latent("SUCCESS", x=0.65, y=0.5, rx=0.08, ry=0.06)
#'
#'  # opwn wall
#'  wall()
#'  title("Drawing three latent variables", col.main="gray20")
#'  
#'  # draw variables
#'  draw(attack)
#'  draw(defense)
#'  draw(success)
#'  }
#'
latent <- 
  function(label = "latent", x = 0.5, y = 0.5, rx = 0.05, ry = 0.05,
           border = "white", lwd = 2, fill = "#5f8bd7", col = "white",
           cex = 1, vfont = NULL, font = 2, family = "sans")
{
  # list of specs
  lat = list(label = label, 
             x = x, 
             y = y, 
             rx = rx, 
             ry = ry,
             north = c(x, y+ry),
             south = c(x, y-ry),
             east = c(x+rx, y), 
             west = c(x-rx, y),
             lcol = border, 
             fill = fill,
             lwd = lwd,
             col = col,
             cex = cex,
             font = font,
             vfont = vfont,
             family = family)
  class(lat) = "latent"
  lat
}
