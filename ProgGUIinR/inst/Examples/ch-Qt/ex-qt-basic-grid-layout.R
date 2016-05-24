
###################################################
### code chunk number 126: LayoutExample
###################################################
window <- Qt$QWidget()
window$setWindowTitle("Layout example")
layout <- Qt$QGridLayout()
window$setLayout(layout)


###################################################
### code chunk number 127: addEntryWidget
###################################################
layout$addWidget(Qt$QLabel("Entry:"), 0, 0)
layout$addWidget(Qt$QLineEdit(), 0, 1, rowspan = 1, colspan=2)


###################################################
### code chunk number 128: addChoiceWidget
###################################################
layout$addWidget(Qt$QLabel("Choice:"), 1, 0)
layout$addWidget(Qt$QComboBox(), 1, 1)


###################################################
### code chunk number 129: addBox
###################################################
layout$addLayout(sub_layout <- Qt$QVBoxLayout(), 
                 1, 2, rowspan=3, 1)
sub_layout$addWidget(label <- Qt$QLabel("Category\nSelector"))
label$setFrameStyle(Qt$QFrame$Box)


###################################################
### code chunk number 130: addLabel
###################################################
layout$addWidget(Qt$QLabel("Text:"), 2, 0, Qt$Qt$AlignTop)
layout$addWidget(edit <- Qt$QTextEdit(), 2, 1)


###################################################
### code chunk number 131: addLabel2
###################################################
layout$addWidget(label <- Qt$QLabel("More info:"), 3, 0, 
                 rowspan = 1, colspan = 2)
label$setSizePolicy(Qt$QSizePolicy$Fixed, 
                    Qt$QSizePolicy$Preferred)
label$setFrameStyle(Qt$QFrame$Box)


###################################################
### code chunk number 132: Layouts.Rnw:433-435
###################################################
layout$setRowStretch(2, 1)                 # third row
layout$setColumnStretch(1,1)               # second column


###################################################
### code chunk number 133: qt-layout-grid-at
###################################################
edit <- layout$itemAtPosition(0, 1)$widget()

##
window$show()
window$raise()
