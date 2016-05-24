# A. Set color and heart parameter values
setwd("C:/aErer")
n <- 500
fun.col <- colorRampPalette(colors = c("red", "white")); fun.col
set.col <- fun.col(n); head(set.col)
set.val <- seq(from = 16, to = 0, length.out = n)

# B. Create a new function to draw one heart as a polygon
oneHeart <- function(r, col) {
  t <- seq(from = 0, to = 2 * pi, length.out = 100)
  x <- r * sin(t) ^ 3
  y <- (13 * r / 16) * cos(t) - (5 * r / 16) * cos(2 * t) - 
       (2 * r / 16) * cos(3 * t) - (r / 16) * cos(4 * t)
  polygon(x, y, col = col, border = NA)
} 

# C. Draw many hearts with mapply()
windows(width = 5.3, height = 3.2); bringToTop(stay = TRUE)
par(mgp = c(0, 0, 0), mai = c(0, 0, 0, 0))
plot.new()
plot.window(xlim = c(-16, 16), ylim = c(-16, 13))
mapply(FUN = oneHeart, set.val, set.col)
showHeart <- recordPlot()

# D. Save the graph on a file device
pdf(file = "fig_showHeart.pdf", width = 5.3, height = 3.2)
replayPlot(showHeart); dev.off()