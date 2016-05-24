
###################################################
### code chunk number 33: Glade.Rnw:16-18
###################################################
builder <- gtkBuilder()
#builder$addFromFile("buildable.xml")
buildable <- system.file("Resources", "buildable.xml", package="ProgGUIinR")
builder$addFromFile(buildable)

###################################################
### code chunk number 34: Glade.Rnw:25-27
###################################################
dialog1 <- builder$getObject("dialog1")
dialog1$showAll()


###################################################
### code chunk number 35: Glade.Rnw:39-43
###################################################
ok_button_clicked <- function(button, userData) {
  message("hello world")
}
builder$connectSignals()

