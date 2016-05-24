###################################################
### code chunk number 182: Dialogs.Rnw:507-508
###################################################
dirname <- ""


###################################################
### code chunk number 183: QFileDialogAPI
###################################################
dialog <- Qt$QFileDialog(NULL, "Choose an R file", getwd(), 
                      name_filter)
dialog$fileMode <- Qt$QFileDialog$ExistingFiles
dialog$setHistory(dirname)


###################################################
### code chunk number 184: Dialogs.Rnw:520-522 (eval = FALSE)
###################################################
if(dialog$exec())
   print(dialog$selectedFiles())

