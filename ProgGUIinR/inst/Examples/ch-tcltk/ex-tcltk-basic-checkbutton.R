
###################################################
### code chunk number 125: setup-window
###################################################
window <- tktoplevel(); tkwm.title(window, "Check button example")
frame <- ttkframe(window); tkpack(frame, expand = TRUE, fill = "both")


###################################################
### code chunk number 126: make-TCL-variables
###################################################
value_var <- tclVar(TRUE)
callback <- function() print(tclvalue(value_var)) # uses global
label_var <- tclVar("check button label")
check_button <- 
  ttkcheckbutton(frame, variable = value_var, 
                 textvariable = label_var, command = callback)
tkpack(check_button)


###################################################
### code chunk number 127: Widgets.Rnw:99-100
###################################################
tkconfigure(check_button, style = "Toolbutton")


###################################################
### code chunk number 128: Widgets.Rnw:112-113
###################################################
tkcget(check_button, "variable" = NULL)

