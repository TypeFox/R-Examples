
###################################################
### code chunk number 151: splitter
###################################################
splitter <- Qt$QSplitter()
splitter$addWidget(Qt$QLabel("One"))
splitter$addWidget(Qt$QLabel("Two"))
splitter$addWidget(Qt$QLabel("Three"))


###################################################
### code chunk number 152: splitterOrientation
###################################################
splitter$setOrientation(Qt$Qt$Vertical)


###################################################
### code chunk number 153: qt-layout-splitter-set
###################################################
splitter$setSizes(c(100L, 200L, 300L))


###################################################
### code chunk number 154: Layouts.Rnw:786-788
###################################################
splitter$show()
splitter$raise()
