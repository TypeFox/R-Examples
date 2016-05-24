## Emilio Torres Manzanera
## University of Oviedo
## Time-stamp: <2014-12-15 15:47 emilio on emilio-despacho>
## ============================================================

##' Creates an XKCD theme
##' 
##' This function creates an XKCD theme
##'
##' @return A layer with the theme. 
##' @import ggplot2
##' @import extrafont
##' @note
##' See the vignette \code{vignette("xkcd-intro")}
##' @export
##' @examples
##' p <- ggplot() + geom_point(aes(mpg, wt), data=mtcars) +
##'      theme_xkcd()
##' p





theme_xkcd <- function(){
  if( "xkcd" %in% extrafont::fonts() ) {
    theme(panel.grid.major = element_blank(),
          ##axis.ticks = element_blank(),
          axis.ticks = element_line(colour = "black"),
          panel.background = element_blank(),
          panel.grid.minor = element_blank(),
          legend.key = element_blank(),
          strip.background = element_blank(),
          text = element_text(size = 16, family = "xkcd"))
  } else {
    warning("Not xkcd fonts installed! See vignette(\"xkcd-intro\")")
    theme(panel.grid.major = element_blank(),
          ##axis.ticks = element_blank(),
          axis.ticks = element_line(colour = "black"),
          panel.background = element_blank(),
          panel.grid.minor = element_blank(),
          legend.key = element_blank(),
          strip.background = element_blank(),
          text = element_text(size = 16))} 
}
