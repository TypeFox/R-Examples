
###################################################
### code chunk number 76: Overview.Rnw:831-845 (eval = FALSE)
###################################################
w <- Qt$QWidget()
w$windowTitle <- "Using Style Sheets"
lyt <- Qt$QHBoxLayout()
w$setLayout(lyt)

b <- Qt$QPushButton("Style sheet example")
lyt$addWidget(b)
b1 <- Qt$QPushButton("Style sheet example")
qsetStyleSheet(color = "red", background = "white", 
               what = "QPushButton", widget = b1)
lyt$addWidget(b1)

w$show()
w$raise()


###################################################
### code chunk number 77: qt-overview-qsetStyleSheet (eval = FALSE)
###################################################
qsetStyleSheet(color = "red", background = "white", 
               what = "QPushButton", widget = button)
