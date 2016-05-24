
###################################################
### code chunk number 196: WidgetsWithModels.Rnw:403-407
###################################################
model <- rGtkDataFrame(mtcars)
view <- gtkTreeView(model)
selection <- view$getSelection()
selection$setMode("single")


###################################################
### code chunk number 197: WidgetsWithModels.Rnw:415-428
###################################################
column <- gtkTreeViewColumn()
view$insertColumnWithAttributes(0, "title", gtkCellRendererText(), text = 0)
## pack in GUI
scrolled_window <- gtkScrolledWindow()
scrolled_window$add(view)
##
window <- gtkWindow(show=FALSE)
window['title'] <- "Multiple selection example"
window$add(scrolled_window)
window$show()
## some selection
selection$selectPath(gtkTreePathNewFromIndices(3)) # set 
#



###################################################
### code chunk number 198: WidgetsWithModels.Rnw:433-435
###################################################
selected <- selection$getSelected()
with(selected, model$getValue(iter, 0)$value)


###################################################
### code chunk number 199: WidgetsWithModels.Rnw:447-455
###################################################
gSignalConnect(selection, "changed", function(selection) {
  selected_rows <- selection$getSelectedRows()
  if(length(selected_rows$retval)) {
    rows <- sapply(selected_rows$retval, 
                   gtkTreePathGetIndices) + 1L
    selected_rows$model[rows, 1]
  }
})
