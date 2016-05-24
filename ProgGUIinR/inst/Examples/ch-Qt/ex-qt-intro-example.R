

###################################################
### code chunk number 2: ex-qtbase.Rnw:7-8
###################################################
library(qtbase)


###################################################
### code chunk number 3: ex-qtbase.Rnw:26-30
###################################################
window <- Qt$QWidget()
label <- Qt$QLabel("Date:")
edit <- Qt$QLineEdit()
button <- Qt$QPushButton("Ok")


###################################################
### code chunk number 4: ex-qtbase.Rnw:49-50
###################################################
window$windowTitle <- "An example"


###################################################
### code chunk number 5: ex-qtbase.Rnw:67-68
###################################################
edit$setInputMask("0000-00-00")


###################################################
### code chunk number 6: ex-qtbase.Rnw:81-86
###################################################
layout <- Qt$QGridLayout()
layout$addWidget(label, row = 0, column = 0, 
                 rowSpan = 1, columnSpan = 1)
layout$addWidget(edit,   0, 1, 1, 1)
layout$addWidget(button, 1, 1, 1, 1)


###################################################
### code chunk number 7: ex-qtbase.Rnw:90-91
###################################################
window$setLayout(layout)


###################################################
### code chunk number 8: ex-qtbase.Rnw:96-97
###################################################
window$show()


###################################################
### code chunk number 9: ex-qtbase.Rnw:107-109
###################################################
handler <- function()  print(edit$text)
qconnect(button, "clicked", handler)


###################################################
### code chunk number 10: ex-qtbase.Rnw:133-137
###################################################
qsetClass("DateValidator", Qt$QValidator, 
          function(parent = NULL) {
            super(parent)
          })


###################################################
### code chunk number 11: ex-qtbase.Rnw:142-150
###################################################
qsetMethod("validate", DateValidator, function(input, pos) {
  if(!grepl("^[0-9]{4}-[0-9]{1,2}-[0-9]{1,2}$", input)) 
    return(Qt$QValidator$Intermediate)
  else  if(is.na(as.Date(input, format="%Y-%m-%d"))) 
    return(Qt$QValidator$Invalid)
  else 
    return(Qt$QValidator$Acceptable)
})


###################################################
### code chunk number 12: ex-qtbase.Rnw:164-166
###################################################
validator <- DateValidator()
edit$setValidator(validator)

