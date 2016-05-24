
###################################################
### code chunk number 237: CompleterExample
###################################################
class_browser <- Qt$QWidget()
layout <- Qt$QFormLayout()
class_browser$setLayout(layout)

layout$addRow("Class name", class_edit <- Qt$QLineEdit())
layout$addRow("Method name", method_edit <- Qt$QLineEdit())


###################################################
### code chunk number 238: Widgets.Rnw:722-724
###################################################
class_completer <- Qt$QCompleter(ls(Qt))
class_edit$setCompleter(class_completer)


###################################################
### code chunk number 239: Widgets.Rnw:730-739
###################################################
qconnect(class_edit, "editingFinished", function() {
  class_name <- class_edit$text
  if(class_name == "") return()
  class_object <- get(class_name, envir = Qt)
  if(!is.null(class_object)) {
    method_completer <- Qt$QCompleter(ls(class_object()))
    method_edit$setCompleter(method_completer)
  }
})


###################################################
### code chunk number 240: Widgets.Rnw:743-748
###################################################
w <- Qt$QWidget()
w$windowTitle <- "Completion example"
w$setLayout(layout)
w$show()
w$raise()

