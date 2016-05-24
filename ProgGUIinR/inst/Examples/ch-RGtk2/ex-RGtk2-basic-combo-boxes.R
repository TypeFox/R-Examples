
###################################################
### code chunk number 125: ComboBoxExample
###################################################
## An example of two comboboxes where 1 updates the other
require(RGtk2)
data(mtcars); library(MASS); data(Cars93) # need some data frames


###################################################
### code chunk number 126: ex-RGtk2-comboboxes.Rnw:11-12
###################################################
library(ProgGUIinR)                     # for avail_dfs, find_vars


###################################################
### code chunk number 127: Widgets
###################################################
window <- gtkWindow(show = FALSE)
window$setTitle("gtkComboBox example")

df_combo <- gtkComboBoxNewText()
var_combo <- gtkComboBoxNewText()


###################################################
### code chunk number 128: Layout
###################################################
vbox <- gtkVBox(); window$add(vbox)
#
vbox1 <- gtkHBox(); vbox$packStart(vbox1)
vbox1$packStart(gtkLabel("Data frames:"))
vbox1$packStart(df_combo)
#
vbox2 <- gtkHBox(); vbox$packStart(vbox2)
vbox2$packStart(gtkLabel("Variable:"))
vbox2$packStart(var_combo)
vbox2$hide()


###################################################
### code chunk number 129: configureComboBoxes
###################################################
sapply(avail_dfs(), df_combo$appendText)
df_combo$setActive(-1)
#
gSignalConnect(df_combo, "changed", function(df_combo, ...) {
  var_combo$getModel()$clear()
  sapply(find_vars(df_combo$getActiveText()),  
         var_combo$appendText)
  vbox2$show()
})


###################################################
### code chunk number 130: ex-RGtk2-comboboxes.Rnw:56-58
###################################################
## show window
window$show()
