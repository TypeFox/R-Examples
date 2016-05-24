

###################################################
### code chunk number 203: WidgetsWithModels.Rnw:504-505
###################################################
require(MASS)


###################################################
### code chunk number 204: basicSort
###################################################
model <- rGtkDataFrame(Cars93)
sorted_model <- gtkTreeModelSortNewWithModel(model)
view <- gtkTreeView(sorted_model)
mapply(view$insertColumnWithAttributes,
       position = -1,
       title = colnames(model),
       cell = list(gtkCellRendererText()),
       text = seq_len(ncol(model)) - 1)
sapply(seq_len(ncol(model)), function(i)
       view$getColumn(i - 1)$setSortColumnId(i - 1))


###################################################
### code chunk number 205: sort-example
###################################################
f <- function(model, iter1, iter2, user.data) {
  types <- c("Compact", "Small", "Sporty", "Midsize", 
             "Large", "Van")
  column <- user.data
  val1 <- model$getValue(iter1, column)$value
  val2 <- model$getValue(iter2, column)$value
  as.integer(match(val1, types) - match(val2, types))
}
sorted_model$setSortFunc(sort.column.id = 3 - 1, sort.func=f, 
                         user.data = 3 - 1)


###################################################
### code chunk number 206: notShown
###################################################
## basic GUI
sw <- gtkScrolledWindow()
sw$add(view)
w <- gtkWindow(show=FALSE)
w['title'] <- "Example of sortable treeview"
w$add(sw)
w$show()

