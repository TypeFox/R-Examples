### R code from vignette source 'ex-RGtk2-menu-popup.Rnw'

###################################################
### code chunk number 1: "menubar-ex"
###################################################
popup <- gtkMenu()                       # top level
popup$append(gtkMenuItem("cut"))
popup$append(gtkMenuItem("copy"))
popup$append(gtkSeparatorMenuItem())
popup$append(gtkMenuItem("paste"))


###################################################
### code chunk number 2: rgtk2-menus-popup-button
###################################################
button <- gtkButton("Click me with right mouse button")
window <- gtkWindow(); window$setTitle("Popup menu example")
window$add(button)


###################################################
### code chunk number 3: ex-RGtk2-menu-popup.Rnw:22-32
###################################################
gSignalConnect(button, "button-press-event",
  f = function(button, event, menu) {
    if(event$getButton() == 3 ||
       (event$getButton() == 1 && # a mac
        event$getState() == GdkModifierType['control-mask'])) 
      gtkMenuPopup(menu, 
                   button = event$getButton(),
                   activate.time = event$getTime())
    return(FALSE)
  }, data = popup)


###################################################
### code chunk number 4: ex-RGtk2-menu-popup.Rnw:43-48
###################################################
sapply(popup$getChildren(), function(child) {
  if(!inherits(child, "GtkSeparatorMenuItem")) # skip these
    gSignalConnect(child, "activate",
           f = function(child, data) message("replace me"))
})


