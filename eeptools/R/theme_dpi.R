#' a ggplot2 theme developed for PDF and PNG for use at the Wisconsin Department of Public Instruction
#' @description This is a custom ggplot2 theme developed for the Wisconsin 
#' Department of Public Instruction
#' @param base_size  numeric, specify the font size as a numeric value, default is 16
#' @param base_family character, specify the font family, this value is optional
#' @details All values are optional
#' @return A theme object which is a list of attributes applied to a ggplot2 object.
#' @source For more information see https://github.com/hadley/ggplot2/wiki/Themes
#' @seealso his uses \code{\link{unit}} from the grid package extensively. 
#' See also \code{\link{theme_bw}} from the ggplot2 package.
#' @author Jared E. Knowles
#' @import ggplot2
#' @export
#'
#' @examples
#' qplot(mpg, wt, data=mtcars) # standard
#' qplot(mpg, wt, data=mtcars) + theme_dpi()
theme_dpi <- function (base_size = 16, base_family = "") {
  theme_grey(base_size = base_size, base_family = base_family) %+replace% 
    theme(axis.title   = element_text(size=rel(0.8), face="bold"),
          axis.title.y = element_text(vjust=0.35, angle=90),
          axis.text  = element_text(size = rel(0.8)),
          axis.ticks = element_line(colour = "black", size=rel(1.5)), 
          legend.key = element_rect(colour = "grey80"), 
          legend.title = element_text(),
          legend.text  = element_text(),
          panel.background = element_rect(fill = "white", colour = NA), 
          panel.border = element_rect(fill = NA, colour = "grey50"), 
          panel.grid.major = element_line(colour = "grey90", size = 0.2), 
          panel.grid.minor = element_line(colour = "grey96", size = 0.5), 
          strip.background = element_rect(fill = "grey90", colour = "grey50"), 
          strip.background = element_rect(fill = "grey90", colour = "grey50"),
          strip.text = element_text(size = rel(0.9), face="bold"),
          strip.text.x = element_text(size = rel(0.9), face="bold"),
          strip.text.y = element_text(size = rel(0.9), face="bold"),
          legend.text = element_text(),
          legend.title = element_text(),
          panel.margin=grid::unit(0.48, "cm")
          )
  }

#' a ggplot2 theme developed for PDF or SVG maps
#' @description This is a custom ggplot2 theme developed for the Wisconsin 
#' Department of Public Instruction for making PDF maps
#' @param base_size numeric, specify the font size, default is 14
#' @param base_family character, specify the font family, this value is optional
#'
#' @return A theme object which is a list of attributes applied to a ggplot2 object.
#' @details All values are optional
#' @source For more information see https://github.com/hadley/ggplot2/wiki/Themes
#' @seealso his uses \code{\link{unit}} from the grid package extensively. 
#' See also \code{\link{theme_bw}} from the ggplot2 package.
#' @author Jared E. Knowles
#' @import ggplot2
#' @export
#'
#' @examples
#' # Data
#' crimes <- data.frame(state = tolower(rownames(USArrests)), USArrests)
#' require(reshape) # for melt
#' crimesm <- melt(crimes, id = 1)
#' # No DPI theme
#' states_map <- map_data("state")
#' ggplot(crimes, aes(map_id = state)) + geom_map(aes(fill = Murder), map = states_map) + 
#'     expand_limits(x = states_map$long, y = states_map$lat)+ labs(title="USA Crime")
#' # Draw map
#' last_plot() + coord_map()
#' # DPI theme
#' ggplot(crimesm, aes(map_id = state)) + geom_map(aes(fill = value), map = states_map) + 
#'     expand_limits(x = states_map$long, y = states_map$lat) + facet_wrap( ~ variable)+theme_dpi_map()
theme_dpi_map <- function(base_size = 14, base_family = "") {
  theme_grey(base_size = base_size, base_family = base_family) %+replace% 
    theme(axis.title = element_blank(),
          axis.text  = element_blank(),
          axis.ticks = element_blank(), 
          legend.key = element_rect(colour = "grey80"), 
          legend.key.size = grid::unit(.8, "lines"),
          legend.title = element_text(size=base_size * 0.8, face="bold"),
          legend.text = element_text(),
          legend.position = "bottom", 
          legend.direction = NULL, 
          legend.justification = "center", 
          panel.background = element_rect(fill = "white",colour = NA), 
          panel.border = element_rect(fill = NA, colour = "grey50"), 
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(), 
          strip.background = element_rect(fill = "grey90", colour = "grey50"), 
          strip.background = element_rect(fill = "grey90", colour = "grey50"),
          strip.text = element_text(size = rel(0.9), face="bold"),
          strip.text.x = element_text(size = rel(0.8), face="bold"),
          strip.text.y = element_text(size = rel(0.8), face="bold"),
          legend.text = element_text(size=base_size * 0.65),
          panel.margin = grid::unit(0.48, "cm"),
          plot.title=element_text(size=base_size * 1.3)
    )
}


#' an alternate ggplot2 theme developed for PDF or SVG maps
#' @description This is a custom ggplot2 theme developed for the Wisconsin 
#' Department of Public Instruction for making PDF maps
#' @param base_size numeric, specify the font size, default is 14
#' @param base_family character, specify the font family, this value is optional
#'
#' @return A theme object which is a list of attributes applied to a ggplot2 object.
#' @details All values are optional
#' @source For more information see https://github.com/hadley/ggplot2/wiki/Themes
#' @seealso his uses \code{\link{unit}} from the grid package extensively. 
#' See also \code{\link{theme_bw}} from the ggplot2 package.
#' @author Jared E. Knowles
#' @import ggplot2
#' @export
#' @examples
#' # Data
#' crimes <- data.frame(state = tolower(rownames(USArrests)), USArrests)
#' require(reshape) # for melt
#' crimesm <- melt(crimes, id = 1)
#' # No DPI theme
#' states_map <- map_data("state")
#' ggplot(crimes, aes(map_id = state)) + geom_map(aes(fill = Murder), map = states_map) + 
#'     expand_limits(x = states_map$long, y = states_map$lat)+ labs(title="USA Crime")
#' # Draw map
#' last_plot() + coord_map()
#' # DPI theme
#' ggplot(crimesm, aes(map_id = state)) + geom_map(aes(fill = value), map = states_map) + 
#'     expand_limits(x = states_map$long, y = states_map$lat) + 
#'     facet_wrap( ~ variable) + theme_dpi_map2()
theme_dpi_map2 <- function(base_size = 14, base_family = "") {
  theme_grey(base_size = base_size, base_family = base_family) %+replace% 
    theme(axis.title = element_blank(),
          axis.text  = element_blank(),
          axis.ticks = element_blank(), 
          legend.background = element_rect(colour="black"),
          legend.key = element_rect(colour = "grey80"), 
          legend.key.size = grid::unit(.7, "lines"),
          legend.title = element_text(size=base_size*0.7, face="bold"),
          legend.text = element_text(size = base_size * 0.5),
          legend.position = c(.13, .15), 
          legend.direction = 'vertical', 
          legend.justification = "center", 
          panel.background = element_rect(fill = "white", colour = NA), 
          panel.border = element_rect(fill = NA, colour = "grey50"), 
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(), 
          strip.background = element_rect(fill = "grey90", colour = "grey50"), 
          strip.background = element_rect(fill = "grey90", colour = "grey50"),
          strip.text = element_text(size = rel(0.9), face="bold"),
          strip.text.x = element_text(size = rel(0.8), face="bold"),
          strip.text.y = element_text(size = rel(0.8), face="bold"),
          legend.text = element_text(size=base_size*0.65),
          panel.margin= grid::unit(0.48, "cm"),
          plot.margin = grid::unit(c(1, 1, 0.5, 0.5), "lines"),
          plot.title=element_text(size = base_size * 1.3)
    )
}


#' an alternate ggplot2 theme developed for PNG or JPG maps
#' @description This is a custom ggplot2 theme developed for the Wisconsin 
#' Department of Public Instruction for making PNG or JPG maps
#' @param base_size numeric, specify the font size, default is 18
#' @param base_family character, specify the font family, this value is optional
#'
#' @return A theme object which is a list of attributes applied to a ggplot2 object.
#' @details All values are optional
#' @source For more information see https://github.com/hadley/ggplot2/wiki/Themes
#' @seealso his uses \code{\link{unit}} from the grid package extensively. 
#' See also \code{\link{theme_bw}} from the ggplot2 package.
#' @author Jared E. Knowles
#' @import ggplot2
#' @export
#' @examples
#' # Data
#' crimes <- data.frame(state = tolower(rownames(USArrests)), USArrests)
#' require(reshape) # for melt
#' crimesm <- melt(crimes, id = 1)
#' # No DPI theme
#' states_map <- map_data("state")
#' ggplot(crimes, aes(map_id = state)) + geom_map(aes(fill = Murder), map = states_map) + 
#'     expand_limits(x = states_map$long, y = states_map$lat)+ labs(title="USA Crime")
#' # Draw map
#' last_plot() + coord_map()
#' # DPI theme
#' ggplot(crimesm, aes(map_id = state)) + geom_map(aes(fill = value), map = states_map) + 
#'     expand_limits(x = states_map$long, y = states_map$lat) + 
#'     facet_wrap( ~ variable) + theme_dpi_mapPNG()
theme_dpi_mapPNG<-function (base_size = 18, base_family = "") {
  theme_grey(base_size = base_size, base_family = base_family) %+replace% 
    theme(axis.title = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(), 
          legend.background = element_rect(colour="black"),
          legend.key = element_rect(colour = "grey80"), 
          legend.key.size = grid::unit(.7, "lines"),
          legend.title = element_text(size = base_size * 1,face="bold"),
          legend.text = element_text(size = base_size * 0.9),
          legend.position = c(.13,.15), 
          legend.margin = grid::unit(0.4, "cm"), 
          legend.direction = 'vertical', 
          legend.justification = "center", 
          panel.background = element_rect(fill = "white", colour = NA), 
          panel.border = element_rect(fill = NA, colour = "grey50"), 
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(), 
          strip.background = element_rect(fill = "grey90", colour = "grey50"), 
          strip.background = element_rect(fill = "grey90", colour = "grey50"),
          strip.text = element_text(size = rel(0.9), face="bold"),
          strip.text.x = element_text(size = rel(0.8), face="bold"),
          strip.text.y = element_text(size = rel(0.8), face="bold"),
          legend.text = element_text(size=base_size * 0.65),
          panel.margin= grid::unit(0.48, "cm"),
          plot.margin = grid::unit(c(1, 1, 0.5, 0.5), "lines"),
          plot.title = element_text(family = base_family, 
                                    size = base_size * 1.4)
    )
}
