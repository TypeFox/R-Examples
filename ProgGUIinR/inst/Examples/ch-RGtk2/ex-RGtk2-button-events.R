

###################################################
### code chunk number 94: CallbackExampleForButton
###################################################
window <- gtkWindow(); button <- gtkButton("click me");
window$add(button)

gSignalConnect(button, "button-press-event", # just mouse
               f = function(widget, event, data) {
                 print(event$getButton())    # which button
                 return(FALSE)               # propagate
               })
gSignalConnect(button, "clicked",            # keyboard too
               f = function(widget, ...) {
                 print("clicked")
               })
