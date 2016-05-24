
###################################################
### code chunk number 76: gtk-container-scrolled-device
###################################################
library(cairoDevice)
device <- gtkDrawingArea()
device$setSizeRequest(600, 400)
asCairoDevice(device)


###################################################
### code chunk number 77: gtk-container-scrolled-construct
###################################################
scrolled <- gtkScrolledWindow()
scrolled$addWithViewport(device)


###################################################
### code chunk number 78: gtk-container-scrolled-zoom
###################################################
zoomPlot <- function(x = 2.0) {
  allocation <- device$getAllocation()$allocation
  device$setSizeRequest(allocation$width * x, 
                        allocation$height * x)
  updateAdjustment <- function(adjustment) {
    adjustment$setValue(x * adjustment$getValue() + 
                        (x - 1) * adjustment$getPageSize()/2)
  }
  updateAdjustment(scrolled$getHadjustment())
  updateAdjustment(scrolled$getVadjustment())
}


###################################################
### code chunk number 79: gtk-container-scrolled-key-press
###################################################
gSignalConnect(scrolled, "key-press-event", 
               function(scrolled, event) {
                 key <- event[["keyval"]]
                 if (key == GDK_plus)
                   zoomPlot(2.0)
                 else if (key == GDK_minus)
                   zoomPlot(0.5)
                 TRUE
               })



###################################################
### code chunk number 80: gtk-container-scrolled-window
###################################################
win <- gtkWindow(show = FALSE)
win$add(scrolled)
win$showAll()


###################################################
### code chunk number 81: gtk-container-scrolled-plot
###################################################
plot(mpg ~ hp, data = mtcars)

