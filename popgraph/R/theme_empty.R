#' A blank theme for plotting networks
#' 
#' This function defines a blank theme for plotting 
#'  graph objects in ggplot.
#' @param base_size The base size of the font
#' @param base_family The font family to use
#' @return A ggplot theme blank and transpatent for plotting
#'  in another program.
#' @author Rodney J. Dyer <rjdyer@@vcu.edu>
#' @export
#' @examples
#' data(lopho)
#' require(ggplot2)
#' require(igraph)
#' layout <- layout.fruchterman.reingold( lopho )
#' V(lopho)$x <- layout[,1]
#' V(lopho)$y <- layout[,2]
#' p <- ggplot() + geom_edgeset( aes(x,y), lopho)
#' p <- p + geom_nodeset( aes(x,y), lopho )
#' p 
#' p + theme_empty()
#' p <- ggplot() + geom_edgeset( aes(x,y,color=weight), lopho)
#' p
#' 
theme_empty <- function (base_size = 12, base_family = "Helvetica"){
    theme(
      line = element_blank(), 
      rect = element_blank(), 
      axis.line = element_blank(), 
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      axis.title = element_blank(), 
      axis.title.y = element_text(angle = 90, vjust = 0.5), 
      axis.ticks.length = grid::unit(0.3, "lines"), 
      axis.ticks.margin = grid::unit(0.5, "lines"),   
      panel.background = element_rect(fill = "transparent", colour = NA), 
      panel.border = element_blank(), 
      panel.grid.major = element_blank(), 
      panel.grid.minor = element_blank(), 
      panel.margin = grid::unit(0.25, "lines"), 
      strip.background = element_blank(), 
      strip.text = element_blank(),
      strip.text.y = element_text(angle = -90), 
      plot.background = element_rect(colour = 'NA', fill = 'transparent'), 
      plot.title = element_text(size = base_size * 1.2), 
      plot.margin = grid::unit(c(1, 1, 0.5, 0.5), "lines"),   
      complete = TRUE)   # denotes that this is a complete theme function
}




