
###################################################
### code chunk number 162: qt-dialogs-input-get-text
###################################################
text <- Qt$QInputDialog$getText(parent = NULL, 
                        title = "Gather text",
                        label = "Enter some text")


###################################################
### code chunk number 163: qt-dialogs-input-get-range
###################################################
even_integer <- Qt$QInputDialog$getInt(parent = NULL, 
                       title="Gather integer",
                       label="Enter an integer from 1 to 10",
                       value=0, min = 2, max = 10, step = 2)


###################################################
### code chunk number 164: qt-dialogs-input-get-item
###################################################
item <- Qt$QInputDialog$getItem(parent = NULL, 
                        title = "Select item",
                        label = "Select a letter",
                        items = LETTERS, current = 17)


###################################################
### code chunk number 165: qt-dialogs-input-explicit
###################################################
dialog <- Qt$QInputDialog()
dialog$setWindowTitle("Select item")
dialog$setLabelText("Select a letter")
dialog$setComboBoxItems(LETTERS)
dialog$setTextValue(LETTERS[18])
dialog$setOptions(Qt$QInputDialog$UseListViewForComboBoxItems)


###################################################
### code chunk number 166: Dialogs.Rnw:245-247
###################################################
if (dialog$exec())
  print(dialog$textValue())
