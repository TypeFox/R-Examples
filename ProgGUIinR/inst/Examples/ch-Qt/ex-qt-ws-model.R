### R code from vignette source 'ex-qt-ws-model.Rnw'

###################################################
### code chunk number 1: ws_model
###################################################
library(qtbase)


###################################################
### code chunk number 2: ex-qt-ws-model.Rnw:16-20
###################################################
qsetClass("WSWatcher", Qt$QObject, function(parent = NULL) {
  super(parent)
  updateVariables()
})


###################################################
### code chunk number 3: ex-qt-ws-model.Rnw:26-27
###################################################
library(digest)


###################################################
### code chunk number 4: ex-qt-ws-model.Rnw:32-33
###################################################
qsetProperty("digests", WSWatcher)


###################################################
### code chunk number 5: objects_changed_property
###################################################
qsetSignal("objectsChanged", WSWatcher)


###################################################
### code chunk number 6: ex-qt-ws-model.Rnw:49-50
###################################################
qsetProperty("objects", WSWatcher, notify = "objectsChanged")


###################################################
### code chunk number 7: ex-qt-ws-model.Rnw:54-56
###################################################
qsetProperty("old_digests", WSWatcher)
qsetProperty("old_objects", WSWatcher)


###################################################
### code chunk number 8: update_variables
###################################################
qsetMethod("updateVariables", WSWatcher, function() {
  x <- sort(ls(envir = .GlobalEnv))
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
### code chunk number 9: change_add
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
### code chunk number 10: addTaskCallback
###################################################
watcher <- WSWatcher()                          # an instance
addTaskCallback(function(expr, value, ok, visible) {
  watcher$updateVariables()
  TRUE
})


###################################################
### code chunk number 11: ex-qt-ws-model.Rnw:111-115 (eval = FALSE)
###################################################
## timer <- Qt$QTimer()
## timer$setSingleShot(FALSE)              # or TRUE for run once
## qconnect(timer,"timeout",function() watcher$updateVariables())
## timer$start(as.integer(3*1000))         # 3 seconds


###################################################
### code chunk number 12: connectSignal
###################################################
qconnect(watcher, "objectsChanged", function() 
         message("workspace objects were updated"))
new_object <- "The change should be announced"


