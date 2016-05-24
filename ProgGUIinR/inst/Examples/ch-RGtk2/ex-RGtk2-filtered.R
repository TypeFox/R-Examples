### R code from vignette source 'ex-RGtk2-filtered.Rnw'

###################################################
### code chunk number 1: ex-RGtk2-filtered.Rnw:5-6
###################################################
library(RGtk2)


###################################################
### code chunk number 2: ex-RGtk2-filtered.Rnw:20-23
###################################################
DF <- data.frame(state.name)
DF$visible <- rep(TRUE, nrow(DF))
model <- rGtkDataFrame(DF)


###################################################
### code chunk number 3: ex-RGtk2-filtered.Rnw:28-31
###################################################
filtered_model <- model$filter()
filtered_model$setVisibleColumn(ncol(DF) - 1)      # offset
view <- gtkTreeView(filtered_model)


###################################################
### code chunk number 4: ex-RGtk2-filtered.Rnw:35-37
###################################################
view$insertColumnWithAttributes(0, "Col", 
                 gtkCellRendererText(), text = 0)


###################################################
### code chunk number 5: ex-RGtk2-filtered.Rnw:45-52
###################################################
entry <- gtkEntry()
gSignalConnect(entry, "changed", function(entry, user.data) {
  pattern <- entry$getText()
  DF <- user.data$getModel()
  values <- DF[, "state.name"]
  DF[, "visible"] <- grepl(pattern, values)
}, data=filtered_model)


###################################################
### code chunk number 6: ex-RGtk2-filtered.Rnw:58-74
###################################################
## not shown, but this places widgets into a simple GUI
window <- gtkWindow(show=FALSE)
window['title'] <- "A filtered data model"
window$setSizeRequest(width=300, height=400)

box <- gtkVBox()
window$add(box)
box$packStart(entry, expand=FALSE)

## add scroll window
sw <- gtkScrolledWindow()
sw$setPolicy("automatic", "automatic")
sw$add(view)
box$packStart(sw, expand=TRUE, fill=TRUE)

window$show()


