
###################################################
### code chunk number 41: Ctrl-q-binding
###################################################
window <- tktoplevel()
button <- ttkbutton(window, text = "Some widget with focus")
tkpack(button)
tkbind(window, "<Control-q>", function() tkdestroy(window))
