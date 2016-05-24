### R code from vignette source 'ex-RGtk2-color-tool-button.Rnw'

###################################################
### code chunk number 1: rgtk2-mennus-toolbar-color-button
###################################################
gdk_color <- gdkColorParse(palette()[1])$color
color_button <- gtkColorButton(gdk_color)


###################################################
### code chunk number 2: rgtk2-menus-toolbar-color-menu
###################################################
colorMenuItem <- function(color) {
  drawing_area <- gtkDrawingArea()
  drawing_area$setSizeRequest(20, 20)
  drawing_area$modifyBg("normal", color)
  image_item <- gtkImageMenuItem(color)
  image_item$setImage(drawing_area)
  image_item
}
color_items <- sapply(palette(), colorMenuItem)
color_menu <- gtkMenu()
for (item in color_items)
  color_menu$append(item)


###################################################
### code chunk number 3: rgtk2-menus-toolbar-color-cb
###################################################
colorMenuItemActivated <- function(item) {
  color <- gdkColorParse(item$getLabel())$color
  color_button$setColor(color)
}
sapply(color_items, gSignalConnect, "activate", 
       colorMenuItemActivated)


###################################################
### code chunk number 4: ex-RGtk2-color-tool-button.Rnw:54-55
###################################################
toolbar <- gtkToolbar()


###################################################
### code chunk number 5: rgtk2-menus-toolbar-menu
###################################################
menu_button <- gtkMenuToolButton(color_button, "Color")
menu_button$setMenu(color_menu)
toolbar$add(menu_button)


## add to GUI
window <- gtkWindow()
window$add(toolbar)
