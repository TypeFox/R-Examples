
###################################################
### code chunk number 85: ex-RGtk2-dialog-layout.Rnw:4-6
###################################################
## layout a basic dialog -- center align
library(RGtk2)


###################################################
### code chunk number 86: gtk-container-table-construct
###################################################
table <- gtkTable(rows = 3, columns = 2, homogeneous = FALSE)


###################################################
### code chunk number 87: ex-RGtk2-dialog-layout.Rnw:26-40
###################################################
size_label <- gtkLabel("Sample size:")
size_combo <- gtkComboBoxNewText()
sapply(c(5, 10, 15, 30), size_combo$appendText)
##
diag_label <- gtkLabel("Diagnostic:")
diag_radio <- gtkVBox()
radiogp <- list()
radiogp$t <- gtkRadioButton(label = "t-statistic")
radiogp$mean <- gtkRadioButton(radiogp, label = "mean")
radiogp$median <- gtkRadioButton(radiogp, label = "median")
sapply(radiogp, diag_radio$packStart)
##
submit_vbox <- gtkVBox()
submit_vbox$packEnd(gtkButton("Run simulation"), expand=FALSE)


###################################################
### code chunk number 88: gtk-container-layout-align
###################################################
size_label['xalign'] <- 1
diag_label['xalign'] <- 1; diag_label['yalign'] <- 0
diag_align <- gtkAlignment(xalign = 0)
diag_align$add(diag_radio)


###################################################
### code chunk number 89: ex-RGtk2-dialog-layout.Rnw:80-94
###################################################
table$attach(size_label, left.attach = 0,1, top.attach = 0,1, 
             xoptions = c("expand", "fill"), yoptions = "")
table$attach(size_combo, left.attach = 1,2, top.attach = 0,1, 
             xoptions = "fill", yoptions = "")
##
table$attach(diag_label, left.attach = 0,1, top.attach = 1,2, 
             xoptions = c("expand", "fill"), 
             yoptions = c("expand", "fill"))
##
table$attach(diag_align, left.attach = 1,2, top.attach = 1,2, 
             xoptions = c("expand", "fill"), yoptions = "")
##
table$attach(submit_vbox, left.attach = 1,2, top.attach = 2,3, 
             xoptions = "", yoptions = c("expand", "fill"))


###################################################
### code chunk number 90: gtk-container-table-spacing
###################################################
table$setColSpacing(0, 5)

