
###################################################
### code chunk number 156: qt-dialogs-construct
###################################################
dialog <- Qt$QMessageBox(icon = Qt$QMessageBox$Warning,
                         title = "Warning!",
                         text = "Warning text...",
                         buttons = Qt$QMessageBox$Ok,
                         parent = NULL)


###################################################
### code chunk number 157: qt-dialogs-extra-text
###################################################
dialog$informativeText <- "Less important warning information"
dialog$detailedText <- "Extra details most do not care to see"


###################################################
### code chunk number 158: qt-dialogs-exec
###################################################
dialog$exec()                       # returns response code
