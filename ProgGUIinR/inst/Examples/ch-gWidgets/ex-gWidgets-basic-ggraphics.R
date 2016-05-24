
###################################################
### code chunk number 140: ggraphicsExample (eval = FALSE)
###################################################
library(gWidgets); options(guiToolkit = "RGtk2")
window <- gwindow("ggraphics example", visible = FALSE)
plot_device <- ggraphics(cont = window)
x <- mtcars$wt; y <- mtcars$mpg
# Identify points
addHandlerClicked(plot_device, handler = function(h, ...) {
  cat(sprintf("You clicked %.2f x %.2f\n", h$x, h$y))
})
# Identify a region, then the points
addHandlerChanged(plot_device, handler = function(h,...) {
  rx <- h$x; ry <- h$y
  if(diff(rx) > diff(range(x))/100 && 
     diff(ry) > diff(range(y))/100) {
    ind <- rx[1] <= x & x <= rx[2] & ry[1] <=y & y <= ry[2]
    if(any(ind))
      print(cbind(x = x[ind], y = y[ind]))
  }
})
visible(window) <- TRUE
#
plot(x, y)
