############################
### code chunk number 139: qt-layout-notebook
###################################################
notebook <- Qt$QTabWidget()


###################################################
### code chunk number 140: qt-layout-notebook-addtab
###################################################
notebook$addTab(Qt$QPushButton("page 1"), "page 1")
icon <- Qt$QIcon("small-R-logo.jpg")
notebook$addTab(Qt$QPushButton("page 2"), icon,  "page 2")


###################################################
### code chunk number 141: qt-layout-notebook-tooltip
###################################################
notebook$setTabToolTip(0, "This is the first page")


###################################################
### code chunk number 142: qt-layout-notebook-current
###################################################
notebook$currentIndex <- 1


###################################################
### code chunk number 143: qt-layout-notebook-pos
###################################################
notebook$tabPosition <- Qt$QTabWidget$South


###################################################
### code chunk number 144: qt-layout-notebook-features
###################################################
notebook$tabsClosable <- TRUE
qconnect(notebook, "tabCloseRequested", function(index) {
  notebook$removeTab(index)
})
notebook$movable <- TRUE
notebook$usesScrollButtons <- TRUE
