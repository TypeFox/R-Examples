### R code from vignette source 'ex-layout.Rnw'

###################################################
### code chunk number 1: ex-layout.Rnw:1-3
###################################################
## layout example
library(qtbase)


###################################################
### code chunk number 2: ex-layout.Rnw:22-26
###################################################
window <- Qt$QWidget()
window$setWindowTitle("Layout example")
grid_layout <- Qt$QGridLayout()
window$setLayout(grid_layout)


###################################################
### code chunk number 3: mainComponents
###################################################
fake_table <- Qt$QLabel("Table widget") 
notebook <- Qt$QTabWidget()
button_layout <- Qt$QHBoxLayout()


###################################################
### code chunk number 4: qt-layout-ex-add
###################################################
grid_layout$addWidget(fake_table, row=0, column=0, 
                     rowspan=1, colspan=1)
grid_layout$addWidget(notebook, 0, 1)
grid_layout$addLayout(button_layout, 1, 1)


###################################################
### code chunk number 5: layoutButtons
###################################################
b <- sapply(c("OK", "Cancel", "Help"), 
            function(i) Qt$QPushButton(i))
button_layout$setDirection(Qt$QBoxLayout$RightToLeft) 
button_layout$addStretch(1L)             # stretch
button_layout$addWidget(b$OK)
button_layout$addWidget(b$Cancel)
button_layout$addSpacing(12L)            # spacing
button_layout$addWidget(b$Help)


###################################################
### code chunk number 6: nbLayout
###################################################
notebook_page <- Qt$QWidget()
notebook$addTab(notebook_page, "Tab label")
notebook$setTabToolTip(0, "A notebook page with a form")


###################################################
### code chunk number 7: nbFormLayout
###################################################
form_layout <- Qt$QFormLayout()
notebook_page$setLayout(form_layout)
l <- sapply(c("name", "rank", "snumber"),  Qt$QLineEdit)
form_layout$addRow("Name", l$name)
form_layout$addRow("Rank", l$rank)
form_layout$addRow("Serial number", l$snumber)


###################################################
### code chunk number 8: ex-layout.Rnw:96-98
###################################################
window$setMinimumSize(width=500, height=400)
window$show()


