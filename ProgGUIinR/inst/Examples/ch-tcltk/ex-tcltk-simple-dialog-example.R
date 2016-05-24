
###################################################
### code chunk number 2: simpleExample
###################################################
library(tcltk)
##
window <- tktoplevel()
tkwm.title(window, "Simple dialog")
##
frame <- ttkframe(window, padding = c(3,3,12,12))
tkpack(frame, expand = TRUE, fill = "both")
##
nested_frame <- ttkframe(frame); tkpack(nested_frame)
##
label <- ttklabel(nested_frame, text = "Enter your name:")
tkpack(label, side = "left")
##
text_var <- tclVar("")
entry <- ttkentry(nested_frame, textvariable = text_var)
tkpack(entry)
##
button_frame <- ttkframe(frame)
tkpack(button_frame, anchor = "ne")
button <- ttkbutton(button_frame, text = "Click")
tkpack(button, side = "right")
##
handler <- function() {
  msg <- sprintf("Hello %s", tclvalue(text_var))
  print(msg)
}
tkconfigure(button, command = handler)
