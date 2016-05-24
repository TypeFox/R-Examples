
###################################################
### code chunk number 38: basic-window-label
###################################################
window <- gtkWindow(show = FALSE)
window$setTitle("Hello World")
label <- gtkLabel("Hello World")
window$add(label)

window$show()


###################################################
### code chunk number 39: gtk-container-window-delete
###################################################
gSignalConnect(window, "delete-event", function(event, ...) {
  dialog <- gtkMessageDialog(parent = window, flags = 0, 
                             type = "question",
                             buttons = "yes-no",
                             "Are you sure you want to quit?")
  out <- dialog$run()
  dialog$destroy()
  out != GtkResponseType["yes"]
})

