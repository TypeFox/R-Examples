#' Add brush to plot
#' This adds a rectangle to a ggplot plot indicating the brush position.
#'
#' @param plot plot object
#' @param x x position of brush
#' @param y y position of brush
#' @param width width of brush
#' @param height height of brush
#' @param just which corner of brush should be determined by x and y position
#' @param fill fill colour for brush.  Use ggplot-alpha for alpha blending.
#' @param col outline colour of brush
#' @author Hadley Wickham \email{h.wickham@@gmail.com}
#' @keywords hplot
#' @export
addbrush <- function(
  plot,
  x, y,
  width = 0.5, height = 0.5,
  just = c("left", "top"),
  fill = "transparent",
  col = "black"
) {
  brush <- data.frame(x = x, y = y, width = width, height = height)
  geom_rect(
    aes_string("x", "y", width = "width", height = "height"),
    data = brush,
    justification = just, fill = fill, colour = col
  )
}

#' Remove hidden points
#' Will remove all hidden points from the plot.
#'
#' @param d ddplot object
#' @author Hadley Wickham \email{h.wickham@@gmail.com}
#' @keywords manip
#' @export
removehiddens <- function(d) {
  d$plots <- lapply(d$plots, function(dd) {
    dd$points <- dd$points[!dd$points$hidden, ]
    dd
  })

  d
}

#' Run All Examples
#' Will run all examples within the package
#' @author Barret Schloerke schloerke@@gmail.com
#' @keywords hplot
#' @export
zeeThemAll <- function() {
  example(ggplot.dd)
  example(ggplot.scatmat)
  example(ggplot.parcoords)
  example(ggplot.timeseries)
}


"%||%" <- function(a, b) {
  if (!is.null(a)) a else b
}
