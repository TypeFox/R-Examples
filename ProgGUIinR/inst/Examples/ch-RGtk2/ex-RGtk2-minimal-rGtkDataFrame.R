### R code from vignette source 'ex-RGtk2-minimal-rGtkDataFrame.Rnw'

###################################################
### code chunk number 1: echo-FALSE
###################################################
library(RGtk2)


###################################################
### code chunk number 2: defineSomeData
###################################################
DF <- data.frame(a=c(1:5), b=c(21.2, rnorm(1), 1, NA, NaN))
store <- rGtkDataFrame(DF)
tv <- gtkTreeView(store)


###################################################
### code chunk number 3: minimalSteps
###################################################
vc <- gtkTreeViewColumn()
QT <- tv$insertColumn(vc, 0)                  # first col. of tree view
vc$setTitle("Column 1")
cr <- gtkCellRendererText()
vc$packStart(cr)
vc$addAttribute(cr, "text", 0)          # first col. of store


###################################################
### code chunk number 4: basicGUI
###################################################
w <- gtkWindow(show=FALSE)
w$setSizeRequest(300,300)
w$setTitle("rGtk2DataFrame example")
sw <- gtkScrolledWindow()
sw$setPolicy("automatic", "automatic")
sw$add(tv)
w$add(sw)
w$showAll()


###################################################
### code chunk number 5: ex-RGtk2-minimal-rGtkDataFrame.Rnw:44-49
###################################################
vc <- gtkTreeViewColumn()
tv$insertColumn(vc, 1)
vc$setTitle("Column 2")
cr <- gtkCellRendererText()
vc$packStart(cr)


###################################################
### code chunk number 6: cellDataFunc
###################################################
cellFunc <- function(vc, cr, model, iter, data) {
  curVal <- model$getValue(iter, data - 1)$value
  if(is.nan(curVal))
    curVal <- "NA"
  else if(is.nan(curVal))
    curVal <- "NaN"
  else
    curVal <- sprintf("%.3f", curVal)
  cr['xalign'] <- 1
  cr['text'] <- curVal
}
QT <- vc$setCellDataFunc(cr, func=cellFunc, 
                         func.data=2) # pass in col. no


