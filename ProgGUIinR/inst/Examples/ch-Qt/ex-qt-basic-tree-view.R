
###################################################
### code chunk number 326: qt-mvc-standard-item-model
###################################################
tree_model <- Qt$QStandardItemModel(rows = 0, columns = 1)


###################################################
### code chunk number 327: qt-mvc-standard-item-set-item
###################################################
by(Cars93, Cars93$Manufacturer, function(DF) {
  tree_model$insertRow(tree_model$rowCount())
  manufacturer <- tree_model$index(tree_model$rowCount()-1L,0)
  tree_model$setData(manufacturer, DF$Manufacturer[1])
  tree_model$insertRows(0, nrow(DF), manufacturer)
  tree_model$insertColumn(0, manufacturer)
  for (i in seq_along(DF$Model)) {
    record <- tree_model$index(i-1L, 0, manufacturer)
    tree_model$setData(record, DF$Model[i])
  }
})


###################################################
### code chunk number 328: Widgets-MVC.Rnw:927-928
###################################################
tree_model <- Qt$QStandardItemModel(rows = 0, columns = 1)


###################################################
### code chunk number 329: qt-mvc-standard-item-rewrite
###################################################
by(Cars93, Cars93$Manufacturer, function(DF) {
  manufacturer <- as.character(DF$Manufacturer[1])
  manufacturer_item <- Qt$QStandardItem(manufacturer)
  tree_model$appendRow(manufacturer_item)
  children <- lapply(as.character(DF$Model), Qt$QStandardItem)
  lapply(children, manufacturer_item$appendRow)
})


###################################################
### code chunk number 330: qt-mvc-tree-view
###################################################
tree_view <- Qt$QTreeView()
tree_view$setModel(tree_model)


###################################################
### code chunk number 331: Widgets-MVC.Rnw:950-963
###################################################
## How to see model three ways
view <- Qt$QTreeView()
view$windowTitle <- "QTreeView"
view$headerHidden <- TRUE
view$setModel(tree_model); view$show(); view$raise()
##
view <- Qt$QTableView()
view$windowTitle <- "QTableView"
view$setModel(tree_model); view$show(); view$raise()
##
view <- Qt$QListView()
view$windowTitle <- "QListView"
view$setModel(tree_model); view$show(); view$raise()


###################################################
### code chunk number 332: qt-mvc-tree-view-header-hidden
###################################################
tree_view$headerHidden <- TRUE

