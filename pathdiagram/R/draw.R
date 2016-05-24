#' @title Draw manifest and latent variables
#' 
#' @description
#' Use this function to draw either manifest or latent variables on a plot.
#' 
#' @param obj An object of either class \code{"manifest"} or \code{"latent"}
#' @author Gaston Sanchez
#' @seealso \code{\link{manifest}}, \code{\link{latent}}
#' @export
#' @examples
#'
#'  \dontrun{
#'  # manifest variables
#'  ingredients = list(
#'    eggs = manifest("eggs", x=0.3, y=0.7, width=0.10, height=0.08),
#'    milk = manifest("milk", x=0.3, y=0.6, width=0.10, height=0.08),
#'    flour = manifest("flour", x=0.3, y=0.5, width=0.10, height=0.08),
#'    sugar = manifest("sugar", x=0.3, y=0.4, width=0.10, height=0.08),
#'    butter = manifest("butter", x=0.3, y=0.3, width=0.10, height=0.08)
#'  )
#'
#'  # latent variables
#'  pancakes = latent("PANCAKES", x=0.6, y=0.6, rx=0.09, ry=0.07)
#'  waffles = latent("WAFFLES", x=0.6, y=0.4, rx=0.09, ry=0.07)
#'
#'  # open wall
#'  wall()
#'  
#'  title("Toy Path Diagram", col.main="gray20")
#'  # draw manifest variables
#'  for (i in 1:length(ingredients)) {
#'     draw(ingredients[[i]])
#'  }
#'  
#'  # draw latent variables
#'  draw(pancakes)
#'  draw(waffles)
#'  # draw arrows
#'  for (i in 1:length(ingredients)) {
#'     arrow(from=ingredients[[i]], to=pancakes, start="east", end="west")
#'     arrow(from=ingredients[[i]], to=waffles, start="east", end="west")
#'  }
#'  }
#'  
draw <- 
function(obj)
{
  if (class(obj) == "manifest")
  {
    rect(obj$xl, obj$yb, obj$xr, obj$yt, col=obj$fill, border=obj$border,
         lwd=obj$lwd)
    text(obj$north[1], obj$east[2], labels=obj$label, col=obj$col,
         cex=obj$cex, vfont=obj$vfont, font=obj$font, family=obj$family)    
  }
  if (class(obj) == "latent")
  {
    plotellipse(rx=obj$rx, ry=obj$ry, mid=c(obj$x, obj$y), 
                lcol=obj$lcol, col=obj$fill, lwd=obj$lwd)
    text(obj$x, obj$y, labels=obj$label, col=obj$col,
         cex=obj$cex, vfont=obj$vfont, font=obj$font, family=obj$family)
  }
}
