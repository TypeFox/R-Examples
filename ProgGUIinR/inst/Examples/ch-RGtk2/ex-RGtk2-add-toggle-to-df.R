### R code from vignette source 'ex-RGtk2-add-toggle-to-df.Rnw'

###################################################
### code chunk number 1: ex-RGtk2-add-toggle-to-df.Rnw:11-13
###################################################
## example showing how to add a toggle button on left of data display
library(RGtk2)


###################################################
### code chunk number 2: FixACRANforSweave
###################################################
repos <- getOption("repos")
repos["CRAN"] <- "http://streaming.stat.iastate.edu/CRAN"
options(repos = repos)


###################################################
### code chunk number 3: getUpgradablePackages
###################################################
old_packages <- 
  old.packages()[,c("Package", "Installed", "ReposVer")]
DF <- as.data.frame(old_packages)


###################################################
### code chunk number 4: ex-RGtk2-add-toggle-to-df.Rnw:33-35
###################################################
doUpdate <- function(old_packages) 
  install.packages(old_packages$Package)


###################################################
### code chunk number 5: ex-RGtk2-add-toggle-to-df.Rnw:47-48
###################################################
model <- rGtkDataFrame(cbind(DF, .toggle=rep(FALSE, nrow(DF))))


###################################################
### code chunk number 6: ex-RGtk2-add-toggle-to-df.Rnw:53-66
###################################################
view <- gtkTreeView()
cell_renderer <- gtkCellRendererToggle() # add toggle
view$insertColumnWithAttributes(0, "", cell_renderer, 
                                active = ncol(DF))
cell_renderer['activatable'] <- TRUE
gSignalConnect(cell_renderer, "toggled", 
               function(cell_renderer, path, user.data) {
                 view <- user.data
                 row <- as.numeric(path) + 1
                 model <- view$getModel()
                 n <- dim(model)[2]
                 model[row, n] <- !model[row, n]
               }, data=view)


###################################################
### code chunk number 7: ex-RGtk2-add-toggle-to-df.Rnw:70-72 (eval = FALSE)
###################################################
mapply(view$insertColumnWithAttributes, -1, colnames(DF), 
        list(gtkCellRendererText()), text = seq_along(DF) -1L)


###################################################
### code chunk number 8: ex-RGtk2-add-toggle-to-df.Rnw:76-77
###################################################
view$setModel(model)


###################################################
### code chunk number 9: ex-RGtk2-add-toggle-to-df.Rnw:86-94
###################################################
button <- gtkButton("Update packages")
gSignalConnect(button, "clicked", function(button, data) {
  view <- data
  model <- view$getModel()
  old_packages <- 
    model[model[, ncol(model)], -ncol(model), drop = FALSE]
  doUpdate(old_packages)
}, data=view)


###################################################
### code chunk number 10: ex-RGtk2-add-toggle-to-df.Rnw:100-112
###################################################
window <- gtkWindow(show = FALSE)
window$setTitle("Installed packages that need upgrading")
window$setSizeRequest(300, 300)

vbox <- gtkVBox(); window$add(vbox)
scrolled_window <- gtkScrolledWindow()
vbox$packStart(scrolled_window, expand = TRUE, fill = TRUE)

scrolled_window$add(view)
scrolled_window$setPolicy("automatic", "automatic")
vbox$packStart(button, expand = FALSE)
window$show()


