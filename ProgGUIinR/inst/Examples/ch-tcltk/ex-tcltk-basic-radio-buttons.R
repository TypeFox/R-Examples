###################################################
### code chunk number 133: radio-button-1
###################################################
window <- tktoplevel(); tkwm.title(window, "Radio example")
frame <- ttkframe(window, padding = c(3,3,12,12)); tkpack(frame, expand = TRUE, fill = "both")


###################################################
### code chunk number 134: radio-button-2
###################################################
values <- c("less", "greater", "two.sided")
var <- tclVar(values[3])                # initial value
callback <- function() print(tclvalue(var))
sapply(values, function(i) {
  radio_button <- ttkradiobutton(frame, variable = var, 
                                 text = i, value = i, 
                                 command = callback)
  tkpack(radio_button, side = "top", anchor = "w")
})


