

###################################################
### code chunk number 255: WidgetsWithModels.Rnw:940-952
###################################################
model <- gtkTreeStore("gchararray")
by(Cars93, Cars93$Manufacturer, function(DF) {
  parent_iter <- model$append()
  model$setValue(parent_iter$iter, column = 0, value = 
                 DF$Manufacturer[1])
  sapply(DF$Model, function(car_model) {
    child_iter <- model$append(parent = parent_iter$iter)
    if (is.null(child_iter$retval)) 
      model$setValue(child_iter$iter, column = 0, 
                     value = car_model)
  })
})


###################################################
### code chunk number 256: WidgetsWithModels.Rnw:957-959
###################################################
iter <- model$getIterFromString("0:0")
model$getValue(iter$iter, column = 0)$value




###################################################
### code chunk number 257: rgtk2-mvc-tree-traverse
###################################################
iter <- model$getIterFirst()
values <- NULL
while(iter$retval) {
  child_iter <- model$iterChildren(iter$iter)
  while(child_iter$retval) {
    values <- c(values, model$get(child_iter$iter, 0)[[1]])
    child_iter$retval <- model$iterNext(child_iter$iter)
  }
  iter$retval <- model$iterNext(iter$iter)
}


##
values

