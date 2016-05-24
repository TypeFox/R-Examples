#############################################
### code chunk number 79: LineEditWithError
###################################################
qsetClass("LineEditWithError", Qt$QLineEdit)


###################################################
### code chunk number 80: setError
###################################################
qsetMethod("setError", LineEditWithError, function(msg) {
  file <- system.file("images/cancel.gif", package="gWidgets")
  qsetStyleSheet("background-image" = sprintf("url(%s)", file),
                 "background-repeat" = "no-repeat",
                 "background-position" = "left top",
                 "padding-left" = "20px",
                 widget = this)
  setToolTip(msg)
})


###################################################
### code chunk number 81: clearError
###################################################
qsetMethod("clearError", LineEditWithError, function() {
  setStyleSheet(NULL)
  setToolTip(NULL)
})


###################################################
### code chunk number 82: testOut
###################################################
edit <- LineEditWithError()
edit$text <- "The quick brown fox..."
edit$setError("Replace with better boilerplate text")
edit$clearError()


###################################################
### code chunk number 83: show_raise
###################################################
w <- Qt$QWidget()
lyt <- Qt$QHBoxLayout()
lyt$addWidget(edit)
w$setLayout(lyt)
w$show()
w$raise()

