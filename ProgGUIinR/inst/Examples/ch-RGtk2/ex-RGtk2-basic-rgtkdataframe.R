
###################################################
### code chunk number 185: WidgetsWithModels.Rnw:84-88
###################################################
data(Cars93, package="MASS")             # mix of classes
model <- rGtkDataFrame(Cars93)
model[1, 4] <- 12
model[1, 4]                              # get value


###################################################
### code chunk number 186: WidgetsWithModels.Rnw:105-106
###################################################
model$setFrame(Cars93[1:5, 1:5])


###################################################
### code chunk number 187: rgtk2-mvc-treeview-construc
###################################################
view <- gtkTreeView(model)


###################################################
### code chunk number 188: rgtk2-mvc-insert-column-hardway
###################################################
column <- gtkTreeViewColumn()
column$setTitle("Manufacturer")
cell_renderer <- gtkCellRendererText()
column$packStart(cell_renderer)
column$addAttribute(cell_renderer, "text", 0)
view$insertColumn(column, 0)


###################################################
### code chunk number 189: rgtk2-mvc-insert-column-easyway
###################################################
view$insertColumnWithAttributes(position = -1, 
                                title = "Model", 
                                cell = gtkCellRendererText(), 
                                text = 2 - 1) # second column


###################################################
### code chunk number 190: rgtk2-mvc-insert-all-columns
###################################################
view <- gtkTreeView(model)
mapply(view$insertColumnWithAttributes,  
       position = -1, 
       title = colnames(model), 
       cell = list(gtkCellRendererText()), 
       text = seq_len(ncol(model)) - 1
       )


###################################################
### code chunk number 191: scrollView
###################################################
window <- gtkWindow()
window$setTitle("Tabular view of data frame")
scrolled_window <- gtkScrolledWindow()
window$add(scrolled_window)
scrolled_window$add(view)

