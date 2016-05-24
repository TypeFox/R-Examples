
###################################################
### code chunk number 201: Widgets.Rnw:240-244
###################################################
window <- Qt$QWidget()
group_box <- Qt$QGroupBox("Cylinders:")
layout <- Qt$QVBoxLayout()
window$setLayout(layout)


###################################################
### code chunk number 202: Widgets.Rnw:248-250
###################################################
btn_group <- Qt$QButtonGroup()
btn_group$exclusive <- FALSE


###################################################
### code chunk number 203: Widgets.Rnw:259-268
###################################################
data(Cars93, package="MASS")
cylinders <- levels(Cars93$Cylinders)
sapply(seq_along(cylinders), function(i) {
  button <- Qt$QCheckBox(sprintf("%s Cylinders", cylinders[i]))
  layout$addWidget(button)
  btn_group$addButton(button, i)
})
sapply(btn_group$buttons(), 
       function(button) button$checked <- TRUE)


###################################################
### code chunk number 204: Widgets.Rnw:280-286
###################################################
checked <- sapply(btn_group$buttons(), function(i) i$checked)
if(any(checked)) {
  checked_cyls <- Cars93$Cylinders %in% cylinders[checked]
  message(sprintf("You've selected %d cases", 
                  sum(checked_cyls)))
}


###################################################
### code chunk number 205: Widgets.Rnw:299-305
###################################################
qconnect(btn_group, "buttonClicked(QAbstractButton*)", 
         function(button) {
           msg <- sprintf("Level '%s': %s", 
                          button$text, button$checked)
           message(msg)
})


###################################################
### code chunk number 206: Widgets.Rnw:307-309
###################################################
window$show()
window$raise()

