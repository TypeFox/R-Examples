### R code from vignette source 'ex-RGtk2-combobox.Rnw'

###################################################
### code chunk number 1: ex-RGtk2-combobox.Rnw:5-7
###################################################
## Example of combo box for colors
library(RGtk2)


###################################################
### code chunk number 2: ex-RGtk2-combobox.Rnw:11-12
###################################################
model <- rGtkDataFrame(palette())


###################################################
### code chunk number 3: comboBox
###################################################
combobox <- gtkComboBox(model)
## color
crc <- gtkCellRendererText()
combobox$packStart(crc, expand=FALSE)                
combobox$addAttribute(crc, "cell-background", 0)
crc$width <- 25
## text
crt <- gtkCellRendererText()
crt['xpad'] <- 5                        # give some space
combobox$packStart(crt)
combobox$addAttribute(crt, "text", 0)


###################################################
### code chunk number 4: ex-RGtk2-combobox.Rnw:36-41 (eval = FALSE)
###################################################
## ## display in a window
win <- gtkWindow(show=FALSE)
win$setTitle("Color test")
win$add(combobox)
win$showAll()


