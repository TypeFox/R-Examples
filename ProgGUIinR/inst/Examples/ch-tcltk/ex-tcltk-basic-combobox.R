
###################################################
### code chunk number 158: Widgets.Rnw:326-329
###################################################
window <- tktoplevel(); tkwm.title(window, "Combo box example")
frame <- ttkframe(window, padding = c(3,3,12,12))
tkpack(frame, expand = TRUE, fill = "both")


###################################################
### code chunk number 159: Widgets.Rnw:332-334
###################################################
values <- state.name
var <- tclVar(values[1])              # initial value


###################################################
### code chunk number 160: Widgets.Rnw:338-344
###################################################
combo_box <- ttkcombobox(frame,
                         values = values,
                         textvariable = var,
                         state = "normal",     # or "readonly"
                         justify = "left")
tkpack(combo_box)


###################################################
### code chunk number 161: Combobox-set-values
###################################################
tkconfigure(combo_box, values = tolower(values))


###################################################
### code chunk number 162: combobox-set-length-1 (eval = FALSE)
###################################################
## tkconfigure(combo_box, values = as.tclObj("New York"))


###################################################
### code chunk number 163: Combobox-set
###################################################
tclvalue(var) <- values[2]              # using tcl variable
tkset(combo_box, values[4])             # by value
tcl(combo_box, "current", 4)            # by index


###################################################
### code chunk number 164: Combobox-get
###################################################
tclvalue(var)                           # TCL variable
tkget(combo_box)                        # get subcommand
tcl(combo_box, "current")               # 0-based index


###################################################
### code chunk number 165: combobox-selection-binding
###################################################
tkbind(combo_box, "<<ComboboxSelected>>", function() {
  print(tclvalue(var))
})


###################################################
### code chunk number 166: Combobox-binding-to-return
###################################################
tkbind(combo_box, "<Return>", function(W) {
  val <- tkget(W)
  cat(as.character(val), "\n")
})

