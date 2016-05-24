#' @name gghcl
#' @keywords colour palette ggplot
#' @author Sven E. Templer
#' @title HTML Colours Like ggplot2
#' @description 
#' Calculate HTML colour code from a palette like ggplot2 uses.
#' @details 
#' See \code{?hcl} for explanation of \code{h}, \code{c} and \code{l}.
#' @param n Numeric value to determine size of palette.
#' @param sub Numeric vector with values within range from \code{1}
#' to \code{n} to subset palette.
#' @param h Hue of the colour. Within range of a circle's degrees.
#' @param c Chroma of the colour.
#' @param l Luminance of the colour. Within range from \code{1} to \code{100}.
#' @param \dots Further arguments passed to function \code{hcl}.
#' @return
#' Returns a character vector containing HTML colour code of the
#' standard \code{ggplot} colour palette.
#' @seealso
#' \link{hcl}
#' @examples
#' #
#' 
#' # Plot some palettes:
#' par(mfrow = c(3,1), mai = c(.1,.1,1,.1))
#' p <- matrix(1:10, 10)
#' image(p, col = gghcl(5), axes = FALSE, main ="gghcl(5)")
#' image(p, col = gghcl(10), axes = FALSE, main = "gghcl(10)")
#' image(p, col = gghcl(10, 1:5), axes = FALSE, main ="gghcl(10, 1:5)")
#' # dev.off() # to reset \code{par}
#' 
#' #

#' @export gghcl
gghcl <- function(n, sub = 1:n, h = c(0, 360) + 15, c = 100, l = 65, ...) {
  h <- seq(h[1], h[2], length = n+1)
  hcl(h, c, l, ...)[sub]  
}
