
###################################################
### code chunk number 207: RadioWithList
###################################################
window <- Qt$QGroupBox("Weight:")
radio_buttons <- 
  list(Qt$QRadioButton("Weight < 3000", w),
       Qt$QRadioButton("3000 <= Weight < 4000", w),
       Qt$QRadioButton("4000 <= Weight", w))


###################################################
### code chunk number 208: qt-widget-radio-layout
###################################################
layout <- Qt$QVBoxLayout()
window$setLayout(layout)
sapply(radio_buttons, layout$addWidget)
radio_buttons[[1]]$setChecked(TRUE)


###################################################
### code chunk number 209: qt-widget-radio-checked
###################################################
radio_buttons[[1]]$checked


###################################################
### code chunk number 210: Widgets.Rnw:348-355
###################################################
sapply(radio_buttons, function(button) {
  qconnect(button, "toggled", function(checked) {
    if(checked) {
      message(sprintf("You checked %s.", button$text))
    }
  })
})


###################################################
### code chunk number 211: qt-widget-radio-group
###################################################
btn_group <- Qt$QButtonGroup()
lapply(radio_buttons, btn_group$addButton)


###################################################
### code chunk number 212: qt-widget-radio-group-checked
###################################################
btn_group$checkedButton()$text


###################################################
### code chunk number 213: Widgets.Rnw:379-381
###################################################
window$show()
window$raise()
