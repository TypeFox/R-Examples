### R code from vignette source 'qtbase.Rnw'

###################################################
### code chunk number 1: setup
###################################################
options(width=72)
library(qtbase)
supported <- 
!length(grep("darwin", R.version$platform)) || 
nzchar(Sys.getenv("SECURITYSESSIONID"))


###################################################
### code chunk number 2: syntax (eval = FALSE)
###################################################
## button <- Qt$QPushButton("Press Me!")
## qconnect(button, "pressed", function() print("Pressed"))
## button$show()


###################################################
### code chunk number 3: syntax-real
###################################################
if (supported) {
button <- Qt$QPushButton("Press Me!")
qconnect(button, "pressed", function() print("Pressed"))
button$show()
}


###################################################
### code chunk number 4: Qt (eval = FALSE)
###################################################
## Qt


###################################################
### code chunk number 5: Qt-real
###################################################
if (supported) {
Qt
}


###################################################
### code chunk number 6: libraries-as-environments (eval = FALSE)
###################################################
## head(ls(Qt))
## Qt$QPushButton


###################################################
### code chunk number 7: libraries-as-environments-real
###################################################
if (supported) {
head(ls(Qt))
Qt$QPushButton
}


###################################################
### code chunk number 8: QWidget (eval = FALSE)
###################################################
## button <- Qt$QPushButton("Press Me!")


###################################################
### code chunk number 9: QWidget-real
###################################################
if (supported) {
button <- Qt$QPushButton("Press Me!")
}


###################################################
### code chunk number 10: tr (eval = FALSE)
###################################################
## Qt$QPushButton$tr("Hello World")


###################################################
### code chunk number 11: tr-real
###################################################
if (supported) {
Qt$QPushButton$tr("Hello World")
}


###################################################
### code chunk number 12: show (eval = FALSE)
###################################################
## button$show()


###################################################
### code chunk number 13: show-real
###################################################
if (supported) {
button$show()
}


###################################################
### code chunk number 14: text (eval = FALSE)
###################################################
## button$text
## button$text <- "PUSH ME!"


###################################################
### code chunk number 15: text-real
###################################################
if (supported) {
button$text
button$text <- "PUSH ME!"
}


###################################################
### code chunk number 16: qconnect (eval = FALSE)
###################################################
## qconnect(button, "pressed", function() print("pushed"))


###################################################
### code chunk number 17: qconnect-real
###################################################
if (supported) {
qconnect(button, "pressed", function() print("pushed"))
}


###################################################
### code chunk number 18: qsetClass (eval = FALSE)
###################################################
## qsetClass("PositiveValidator", Qt$QValidator)


###################################################
### code chunk number 19: qsetClass-real
###################################################
if (supported) {
qsetClass("PositiveValidator", Qt$QValidator)
}


###################################################
### code chunk number 20: list-validator-class (eval = FALSE)
###################################################
## PositiveValidator


###################################################
### code chunk number 21: list-validator-class-real
###################################################
if (supported) {
PositiveValidator
}


###################################################
### code chunk number 22: validate (eval = FALSE)
###################################################
## validatePositive <- function(input, pos) {
##   val <- suppressWarnings(as.integer(input))
##   if (!is.na(val)) {
##     if (val > 0)
##     Qt$QValidator$Acceptable
##     else Qt$QValidator$Invalid
##   } else {
##     if (input == "")
##     Qt$QValidator$Acceptable
##     else Qt$QValidator$Invalid
##   }
## }


###################################################
### code chunk number 23: validate-real
###################################################
if (supported) {
validatePositive <- function(input, pos) {
  val <- suppressWarnings(as.integer(input))
  if (!is.na(val)) {
    if (val > 0)
    Qt$QValidator$Acceptable
    else Qt$QValidator$Invalid
  } else {
    if (input == "")
    Qt$QValidator$Acceptable
    else Qt$QValidator$Invalid
  }
}
}


###################################################
### code chunk number 24: qsetMethod (eval = FALSE)
###################################################
## qsetMethod("validate", PositiveValidator, validatePositive)


###################################################
### code chunk number 25: qsetMethod-real
###################################################
if (supported) {
qsetMethod("validate", PositiveValidator, validatePositive)
}


###################################################
### code chunk number 26: construct-validator (eval = FALSE)
###################################################
## validator <- PositiveValidator()


###################################################
### code chunk number 27: construct-validator-real
###################################################
if (supported) {
validator <- PositiveValidator()
}


###################################################
### code chunk number 28: text-entry (eval = FALSE)
###################################################
## e <- Qt$QLineEdit()
## v <- PositiveValidator(e)
## e$setValidator(v)
## e$show()


###################################################
### code chunk number 29: text-entry-real
###################################################
if (supported) {
e <- Qt$QLineEdit()
v <- PositiveValidator(e)
e$setValidator(v)
e$show()
}


###################################################
### code chunk number 30: extend-window-title (eval = FALSE)
###################################################
## qsetClass("SaveConfirmationDialog", Qt$QMessageBox, 
## function(filename = NULL, parent = NULL) 
## {
##   super(icon = Qt$QMessageBox$Question, title = "Save confirmation", 
##   text = "Save the current document?", 
##   buttons = Qt$QMessageBox$Cancel | Qt$QMessageBox$Discard | 
##   Qt$QMessageBox$Save,
##   parent = parent)
##   this$filename <- filename
## })


###################################################
### code chunk number 31: extend-window-title-real
###################################################
if (supported) {
qsetClass("SaveConfirmationDialog", Qt$QMessageBox, 
function(filename = NULL, parent = NULL) 
{
  super(icon = Qt$QMessageBox$Question, title = "Save confirmation", 
  text = "Save the current document?", 
  buttons = Qt$QMessageBox$Cancel | Qt$QMessageBox$Discard | 
  Qt$QMessageBox$Save,
  parent = parent)
  this$filename <- filename
})
}


###################################################
### code chunk number 32: accept-override (eval = FALSE)
###################################################
## qsetMethod("accept", SaveConfirmationDialog, function() {
##   saveDocument(filename)
##   super("accept")
## })


###################################################
### code chunk number 33: accept-override-real
###################################################
if (supported) {
qsetMethod("accept", SaveConfirmationDialog, function() {
  saveDocument(filename)
  super("accept")
})
}


