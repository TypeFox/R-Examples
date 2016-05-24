
###################################################
### code chunk number 161: QMEssageBoxAPI
###################################################
dialog <- Qt$QMessageBox()
dialog$windowTitle <- "[This space for rent]"
dialog$text <- "This is the main text"
dialog$informativeText <- "This should give extra info"
dialog$detailedText <- "And this provides\neven more detail"

dialog$icon <- Qt$QMessageBox$Critical
dialog$standardButtons <- 
  Qt$QMessageBox$Cancel | Qt$QMessageBox$Ok
## 'Cancel' instead of 'Ok' is the default
dialog$setDefaultButton(Qt$QMessageBox$Cancel)
##
if(dialog$exec() == Qt$QMessageBox$Ok) 
  print("A Ok")
