### R code from vignette source 'ex-qt-treewidget.Rnw'

###################################################
### code chunk number 1: ex-qt-treewidget.Rnw:9-11
###################################################
## use QTreeView to make workspace browser
require(qtbase)


###################################################
### code chunk number 2: ex-qt-treewidget.Rnw:24-114
###################################################
## From an earlier example
###################################################
### chunk number 1: ws_watcher
###################################################
library(qtbase)


###################################################
### chunk number 2: 
###################################################
qsetClass("WSWatcher", Qt$QObject, function(parent=NULL) {
  super(parent)
  updateVariables()
})


###################################################
### chunk number 3: 
###################################################
library(digest)


###################################################
### chunk number 4: 
###################################################
qsetProperty("digests", WSWatcher)


###################################################
### chunk number 5: objects_changed_property
###################################################
qsetSignal("objectsChanged", WSWatcher)


###################################################
### chunk number 6: 
###################################################
qsetProperty("objects", WSWatcher, notify="objectsChanged")


###################################################
### chunk number 7: 
###################################################
qsetProperty("old_digests", WSWatcher)
qsetProperty("old_objects", WSWatcher)


###################################################
### chunk number 8: update_variables
###################################################
qsetMethod("updateVariables", WSWatcher, function() {
  x <- sort(ls(envir=.GlobalEnv))
  objs <- sapply(mget(x, .GlobalEnv), digest)

  if((length(objs) != length(digests)) ||
     length(digests) == 0 ||
     any(objs != digests)) {
    this$old_digests <- digests         # old
    this$old_objects <- objects
    this$digests <- objs                # update cache
    this$objects <- x                   # emits signal         
  }
  invisible()
})


###################################################
### chunk number 9: change_add
###################################################
qsetMethod("changedVariables", WSWatcher, function() {
  changed <- setdiff(old_digests, digests)
  old_objects[old_digests %in% changed]
})
##
qsetMethod("addedVariables", WSWatcher, function() {
  added <- setdiff(digests, old_digests)
  objects[digests %in% added]
})





###################################################
### code chunk number 3: ex-qt-treewidget.Rnw:121-141
###################################################
addItem <- function(varname, parent_object, parent_item) {
          
  obj <- parent_object[[varname]]
  ## main interaction with tree model
  item <- Qt$QStandardItem(varname)
  class_item <- Qt$QStandardItem(paste(class(obj), 
                                      collapse = ", "))
  parent_item$appendRow(list(item, class_item))

  ## Recursively create ancestor items, if needed
  nms <- NULL
  if (is.recursive(obj)) {
    if (is.environment(obj))
      nms <- ls(obj)
    else if (!is.null(names(obj)))
      nms <- names(obj)
  }
  sapply(nms, addItem, parent_item = item, 
         parent_object = obj)
}


###################################################
### code chunk number 4: updateTopLevelItems
###################################################
updateTopLevelItems <- function(ws_watcher, view, 
                                env = .GlobalEnv) {
  ## remove these (by index)
  remove <- ws_watcher$changedVariables()
  cur_shown <- sapply(seq(model$rowCount()), 
                 function(i) model$index(i - 1, 0)$data())
  indices_to_remove <- which(cur_shown == remove)
  indices_to_remove <- sort(indices_to_remove, decreasing=TRUE)  
  ## add these (by variable name)
  new_names <- ws_watcher$addedVariables()
  
  ## replace/add these
  model <- view$model()
  view$updatesEnabled <- FALSE
  if(length(indices_to_remove))
    sapply(indices_to_remove -1L, model$removeRow)
  ## add
  sapply(new_names, addItem, parent_object = env,
         parent_item = model$invisibleRootItem())
  model$sort(0, Qt$Qt$AscendingOrder)
  view$updatesEnabled <- TRUE
}


###################################################
### code chunk number 5: ex-qt-treewidget.Rnw:181-192
###################################################
initializeTopLevelItems <- function(ws_watcher, view, 
                                    env = .GlobalEnv) 
{
   current_names <- ws_watcher$objects
   model <- view$model()
   view$updatesEnabled <- FALSE
   sapply(current_names, addItem, parent_object = env, # add
          parent_item = model$invisibleRootItem())
   model$sort(0, Qt$Qt$AscendingOrder)
   view$updatesEnabled <- TRUE
}


###################################################
### code chunk number 6: showTree
###################################################
model <- Qt$QStandardItemModel(rows = 0, columns = 2)
model$setHorizontalHeaderLabels(c("Name", "Class"))
view <- Qt$QTreeView()
view$windowTitle <- "Workspace Browser"
view$headerHidden <- FALSE
view$setModel(model)


###################################################
### code chunk number 7: initialize (eval = FALSE)
###################################################
ws_watcher <- WSWatcher()
ws_watcher$updateVariables()
initializeTopLevelItems(ws_watcher, view)


###################################################
### code chunk number 8: objectsChanged (eval = FALSE)
###################################################
qconnect(ws_watcher, "objectsChanged", function() 
         updateTopLevelItems(ws_watcher, view))


###################################################
### code chunk number 9: taskCallback (eval = FALSE)
###################################################
## ## add callback
id <- addTaskCallback(function(expr, value, ok, visible) {
  ws_watcher$updateVariables()
  TRUE
})
## view
view$show()
view$raise()


