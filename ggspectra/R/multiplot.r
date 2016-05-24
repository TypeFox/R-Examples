#' Multiple plot function
#'
#' Grid based; allows multiple plots arraged in a matrix and \code{print}ed to
#' any R device. ggplot objects can be passed in ..., or to plotlist (as a list
#' of ggplot objects)
#'
#' @param ... one or more ggplot objects
#' @param plotlist list of ggplot objects
#' @param  cols numerical   Number of columns in layout
#' @param layout  A numeric matrix specifying the layout. If present, 'cols' is
#'   ignored.
#'
#' @details ggplot objects can be passed in ..., or to plotlist (as a list of
#'   ggplot objects) If the layout is something like matrix(c(1,2,3,3), nrow=2,
#'   byrow=TRUE), then plot 1 will go in the upper left, 2 will go in the upper
#'   right, and 3 will go all the way across the bottom.
#'
#' @references from http://www.cookbook-r.com/
#'
#' @note Modified from example by Winston Chang found in the Cookbook for R
#' Licenced under CC BY-SA
#'
#' @examples
#' library(photobiology)
#' multiplot(plot(sun.spct), plot(yellow_gel.spct), cols = 2)
#'
#' @export
#'
multiplot <- function(..., plotlist=NULL, cols=1, layout=NULL) {

  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)

  numPlots = length(plots)

  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }

  if (numPlots == 1) {
    print(plots[[1]])

  } else {
    # Set up the page
    grid::grid.newpage()
    grid::pushViewport(grid::viewport(layout = grid::grid.layout(nrow(layout), ncol(layout))))

    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))

      print(plots[[i]], vp = grid::viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}
