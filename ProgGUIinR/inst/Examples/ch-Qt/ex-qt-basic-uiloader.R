## This only works if the library is present

###################################################
### code chunk number 84: has_loader
###################################################
has_uiloader <- !is.null(Qt$QUiLoader)


###################################################
### code chunk number 85: qt-overview-designer-load (eval = FALSE)
###################################################
## loader <- Qt$QUiLoader()
## widget <- loader$load(Qt$QFile("textfinder.ui"))
textfinder.ui <- system.file("Resources", "textfinder.ui", package="ProgGUIinR")

###################################################
### code chunk number 86: Overview.Rnw:964-967
###################################################
if(has_uiloader) {
loader <- Qt$QUiLoader()
widget <- loader$load(Qt$QFile("textfinder.ui"))
}


###################################################
### code chunk number 88: Overview.Rnw:988-991
###################################################
if(has_uiloader) {
find_button <- qfindChild(widget, "findButton") # by name
line_edit <- qfindChild(widget, "lineEdit")
}



###################################################
### code chunk number 90: Overview.Rnw:1001-1004
###################################################
if(has_uiloader) {
qconnect(find_button, "clicked", function() {
  findText(line_edit$text)
})
}



###################################################
### code chunk number 92: Overview.Rnw:1018-1021
###################################################
if(has_uiloader) {
qsetClass("MyMainWindow", Qt$QWidget, function() {
  loader <- Qt$QUiLoader()
  widget <- loader$load(Qt$QFile("textfinder.ui"), this)
  Qt$QMetaObject$connectSlotsByName(this)
})
}


###################################################
### code chunk number 94: Overview.Rnw:1038-1041
###################################################
if(has_uiloader) {
qsetSlot("on_findButton_clicked", MyMainWindow, function() {
  findText(line_edit$text)
})
}



###################################################
### code chunk number 96: Overview.Rnw:1050-1053
###################################################
if(has_uiloader) {
MyMainWindow()
}

