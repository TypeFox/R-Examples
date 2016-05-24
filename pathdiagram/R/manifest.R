#' @title Set specifications of a manifest variable
#' 
#' @description
#' Use this function to specify the graphic characteristics of a manifest 
#' variable. The specifications will be used by the function \code{draw} 
#' to plot manifest variables (in a path diagram).
#' 
#' @details
#' Manifest variables are drawn as rectangles.
#' 
#' @param label A character string with the label to be displayed.
#' @param x x-axis coordinate for center of rectanlge.
#' @param y y-axis coordinate for center of rectangle.
#' @param width width of the rectangle.
#' @param height height of the rectangle.
#' @param border color of the border.
#' @param fill color to fill the rectangle.
#' @param lwd width of the border.
#' @param col color of the label.
#' @param cex numeric character expansion of the label. 
#' @param vfont font family of the label.
#' @param font An integer specifying which font to use for the label.
#' See \code{\link{par}}
#' @param family The name of a font family for drawing text. 
#' Standard values are \code{"serif"}, \code{"sans"} and \code{"mono"}.
#' @return An object of class \code{"manifest"}, which is a list with the
#' specified parameters to draw manifest variables.
#' @author Gaston Sanchez
#' @seealso \code{\link{latent}}, \code{\link{draw}}
#' @export
#' @examples
#'
#'  \dontrun{
#'  # manifest variables
#'  eggs = manifest("eggs", x=0.3, y=0.7, width=NULL, height=0.08)
#'  milk = manifest("milk", x=0.4, y=0.6, width=NULL, height=0.08)
#'  flour = manifest("flour", x=0.5, y=0.5, width=NULL, height=0.08)
#'  sugar = manifest("sugar", x=0.6, y=0.4, width=NULL, height=0.08)
#'  butter = manifest("butter", x=0.7, y=0.3, width=NULL, height=0.08)
#'
#'  # open wall
#'  wall()
#'  title("Five manifest variables", col.main="gray20")
#'  
#'  # draw manifest variables
#'  draw(eggs)
#'  draw(milk)
#'  draw(flour)
#'  draw(sugar)
#'  draw(butter)  
#'  }
#'
manifest <- 
  function(label = "manifest", x = 0.5, y = 0.5, width = NULL, height = 0.1,
           border = "white", fill = "#9dbafa", lwd = 1, col = "gray20", cex = 1,
           vfont = NULL, font = 1, family = "sans")
{
  if (is.null(width)) width = nchar(label) / 50
  xl = x - width / 2
  yb = y - height / 2
  xr = x + width / 2
  yt = y + height / 2
  
  # list of specs
  manif = list(label = label, 
             xl = xl, 
             yb = yb, 
             xr = xr, 
             yt = yt,
             north = c(x, yt),
             south = c(x, yb),
             east = c(xr, y), 
             west = c(xl, y),
             border = border,
             fill = fill,
             lwd = lwd,
             col = col,
             cex = cex,
             vfont = vfont,
             font = font,
             family = family)
  class(manif) = "manifest"
  manif
}
