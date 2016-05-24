### R code from vignette source 'ex-RGtk2-ImageForGraphics.Rnw'

###################################################
### code chunk number 1: ex-RGtk2-ImageForGraphics.Rnw:1-2
###################################################
library(RGtk2)


###################################################
### code chunk number 2: ex-RGtk2-ImageForGraphics.Rnw:12-17
###################################################
window <- gtkWindow(show = FALSE)
window$setTitle("Graphic window")
window$setSizeRequest(400, 400)
hbox <- gtkHBox(); window$add(hbox)
window$showAll()


###################################################
### code chunk number 3: ex-RGtk2-ImageForGraphics.Rnw:25-27
###################################################
theSize <- hbox$getAllocation()$allocation
width <- theSize$width; height <- theSize$height


###################################################
### code chunk number 4: ex-RGtk2-ImageForGraphics.Rnw:33-38
###################################################
require(cairoDevice)
pixmap <- gdkPixmap(drawable = NULL, 
                    width = width, height = height, depth=24)
asCairoDevice(pixmap)
hist(rnorm(100))


###################################################
### code chunk number 5: ex-RGtk2-ImageForGraphics.Rnw:43-45
###################################################
image <- gtkImage(pixmap = pixmap)
hbox$packStart(image, expand = TRUE, fill = TRUE)


