## Example using tkwait to make a dialog modal

###################################################
### code chunk number 118: tkwait
###################################################
msg <- "We care ..."
dialog <- tktoplevel(); tkwm.withdraw(dialog)
tkwm.overrideredirect(dialog, TRUE)   # no decoration
frame <- ttkframe(dialog, padding = 5)
tkpack(frame, expand = TRUE, fill = "both")
tkpack(ttklabel(frame, text = msg), pady = 5)


###################################################
### code chunk number 119: waitVariable
###################################################
flag <- tclVar("")
tkpack(ttkbutton(frame, text="dismiss", command=function() {
  tkgrab.release(dialog)
  tclvalue(flag) <- "Destroy"
}))


###################################################
### code chunk number 120: Dialogs.Rnw:72-74 (eval = FALSE)
###################################################
tkwm.deiconify(dialog)
tkwait.variable(flag)


## The following code won't execute until *after* the dismiss button is clicked causing
## the "flag" variable to be changed.


###################################################
### code chunk number 121: tkdestroy
###################################################
print("Thanks")
tkdestroy(dialog)

