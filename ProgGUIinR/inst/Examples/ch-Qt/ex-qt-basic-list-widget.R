
###################################################
### code chunk number 371: qt-mvc-listwidget-additems
###################################################
list_widget <- Qt$QListWidget()
list_widget$addItems(state.name)

list_widget$show()
list_widget$raise()


###################################################
### code chunk number 372: qt-mvc-listwidget-item
###################################################
item <- Qt$QListWidgetItem("Puerto Rico", list_widget)


###################################################
### code chunk number 373: qt-mvc-listwidget-itemat
###################################################
first <- list_widget$item(0)
first$text()


###################################################
### code chunk number 374: qt-mvc-listwidget-selectionmode
###################################################
list_widget$selectionMode <- Qt$QListWidget$ExtendedSelection


###################################################
### code chunk number 375: qt-mvc-listwidget-select
###################################################
sapply(grep("^A", state.name), 
       function(i) list_widget$item(i - 1)$setSelected(TRUE))


###################################################
### code chunk number 376: qt-mvc-listwidget-selected
###################################################
selected_items <- list_widget$selectedItems()
sapply(selected_items, qinvoke, "text")


###################################################
### code chunk number 377: qt-mvc-listwidget-selectionchanged
###################################################
qconnect(list_widget, "itemSelectionChanged", function() {
  selected <- list_widget$selectedItems()
  selected_text <- sapply(selected, qinvoke, "text")
  message("Selected: ", paste(selected_text, collapse = ", "))
})


###################################################
### code chunk number 378: qt-mvc-listwidget-checked
###################################################
items <- sapply(seq(list_widget$count) - 1L, list_widget$item)
sapply(items, qinvoke, "setCheckState", Qt$Qt$Unchecked)
## check selected
selected <- list_widget$selectedItems()
sapply(selected, function(x) x$setCheckState(Qt$Qt$Checked))
## clear selection now
list_widget$selectionModel()$clear()
list_widget$selectionMode <- Qt$QListWidget$NoSelection


###################################################
### code chunk number 379: Widgets-MVC.Rnw:1225-1227
###################################################
state <- sapply(items, "qinvoke", "checkState")
head(state, n = 8)                     # 2 is checked, 0 not
