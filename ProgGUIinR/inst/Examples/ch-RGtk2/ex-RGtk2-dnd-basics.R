
###################################################
### code chunk number 179: BasicComponents.Rnw:861-868
###################################################
TARGET.TYPE.TEXT   <- 80                # our enumeration
TARGET.TYPE.PIXMAP <- 81                  
widgetTargetTypes <- 
  list(text = gtkTargetEntry("text/plain", 0, 
         TARGET.TYPE.TEXT),
       pixmap = gtkTargetEntry("image/x-pixmap", 0, 
         TARGET.TYPE.PIXMAP))


###################################################
### code chunk number 180: BasicComponents.Rnw:878-886
###################################################
window <- gtkWindow(); window['title'] <- "Drag Source"
drag_source_widget <-  gtkButton("Text to drag")
window$add(drag_source_widget)

gtkDragSourceSet(drag_source_widget,
       start.button.mask=c("button1-mask", "button3-mask"),
       targets=widgetTargetTypes[["text"]],
       actions="copy")


###################################################
### code chunk number 181: BasicComponents.Rnw:899-903
###################################################
gSignalConnect(drag_source_widget, "drag-data-get", 
               function(widget, context, sel, tType, eTime) {
                 sel$setText(widget$getLabel()) 
               })


###################################################
### code chunk number 182: BasicComponents.Rnw:914-922
###################################################
window <- gtkWindow(); window['title'] <- "Drop Target"
drop_target_widget <- gtkButton("Drop here")
window$add(drop_target_widget)

gtkDragDestSet(drop_target_widget,
               flags="all", 
               targets=widgetTargetTypes[["text"]],
               actions="copy")


###################################################
### code chunk number 183: BasicComponents.Rnw:940-945
###################################################
gSignalConnect(drop_target_widget, "drag-data-received", 
       function(widget, context, x, y, sel, tType, eTime) {
         dropdata <- sel$getText()
         widget$setLabel(rawToChar(dropdata))
       })

