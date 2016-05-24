### R code from vignette source 'ex-RGtk2-simple-tree.Rnw'

###################################################
### code chunk number 1: notShown
###################################################
## define tstore, but aslo in earlier example so not shown
data(Cars93, package="MASS")
model <- gtkTreeStore("gchararray")
Manufacturers <- Cars93$Manufacturer
Makes <- split(Cars93[,"Model"], Manufacturers)
for(i in unique(Manufacturers)) {
  piter <- model$append()              # parent
  model$setValue(piter$iter, column=0, value=i)
  for(j in Makes[[i]]) { 
    sibiter <- model$append(parent=piter$iter) # child
    if(is.null(sibiter$retval)) 
      model$setValue(sibiter$iter,column=0, value=j)
  }
}


###################################################
### code chunk number 2: makeView
###################################################
view <- gtkTreeView()
view$insertColumnWithAttributes(0, "Make", 
           gtkCellRendererText(), text = 0)
view$setModel(model)


###################################################
### code chunk number 3: makeGUI
###################################################
w <- gtkWindow(show=FALSE)
w['title'] <- "Example of changing models"
sw <- gtkScrolledWindow()
sw$add(view)
w$add(sw)
w$show()


###################################################
### code chunk number 4: ex-RGtk2-simple-tree.Rnw:45-47
###################################################
model <- rGtkDataFrame(Cars93[,"Model", drop=FALSE])
view$setModel(model)


