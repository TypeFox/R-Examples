
###################################################
### code chunk number 2: Introduction.Rnw:2-3
###################################################
require(RGtk2)


###################################################
### code chunk number 3: gtk-overview-initial-example
###################################################
button <- gtkButton("Click Me")
button['image'] <- gtkImage(stock = "gtk-apply", 
                            size = "button")
gSignalConnect(button, "clicked", function(button) {
  message("Hello World!")
})
##
window <- gtkWindow(show = FALSE)
window$add(button)
window$showAll()

