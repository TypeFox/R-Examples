#' @importFrom graphics lines
lines.barplot <- function(bound.lower, bound.upper, value, bgcolor = "#eeeeee", fgcolor = "#666666", lwd = 8){
  
  stopifnot(is.numeric(bound.lower))
  stopifnot(length(bound.lower) == 2)
  
  stopifnot(is.numeric(bound.upper))
  stopifnot(length(bound.upper) == 2)
  
  stopifnot(is.numeric(value))
  stopifnot(value >= 0 && value <= 1)
  
  x = c(bound.lower[1], bound.upper[1])
  y = c(bound.lower[2], bound.upper[2])
  lines(x, y, col = bgcolor, lwd = lwd, lend = "butt")
  
  value.x <- c(bound.lower[1], bound.lower[1] + value * (bound.upper[1] - bound.lower[1]))
  value.y <- c(bound.lower[2], bound.lower[2] + value * (bound.upper[2] - bound.lower[2]))
  lines(value.x, value.y, col = fgcolor, lwd = lwd, lend = "butt")
}

# plot(0,0)
# lines.barplot(value = 0.5, bound.lower = c(-0.5,0.5), bound.upper = c(0.5,0.5))
# plot(0,0)
# lines.barplot(value = 0.5, bound.lower = c(0.6,-0.5), bound.upper = c(0.6,0.5))
# plot(0,0)
# lines.barplot(value = 0.5, bound.lower = c(-0.6,-0.5), bound.upper = c(+0.6,0.5))
