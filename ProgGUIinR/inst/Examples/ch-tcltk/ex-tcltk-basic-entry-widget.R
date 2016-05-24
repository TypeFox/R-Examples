
###################################################
### code chunk number 135: entryExample
###################################################
window <- tktoplevel()
tkwm.title(window, "Entry widget test")
frame <- ttkframe(window, padding = c(3,3,12,12)); tkpack(frame, expand = TRUE, fill = "both")


###################################################
### code chunk number 136: entryExampleDef
###################################################
txt_var <- tclVar("initial value")
entry <- ttkentry(window, textvariable = txt_var)
tkpack(entry)


###################################################
### code chunk number 137: Widgets.Rnw:227-229
###################################################
tclvalue(txt_var)
tclvalue(txt_var) <- "set value"


###################################################
### code chunk number 138: tkget
###################################################
tkget(entry)


###################################################
### code chunk number 139: tkinsert
###################################################
tkinsert(entry, "end", "new text")


###################################################
### code chunk number 140: Widgets.Rnw:253-254
###################################################
tkdelete(entry, 0, 4)


###################################################
### code chunk number 141: Widgets.Rnw:262-263
###################################################
tkicursor(entry, 0)                         # move to beginning


###################################################
### code chunk number 142: Widgets.Rnw:271-272
###################################################
tkselection.range(entry, 0, "end")

