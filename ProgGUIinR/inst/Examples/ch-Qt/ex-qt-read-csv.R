### R code from vignette source 'ex-qt-read-csv.Rnw'

###################################################
### code chunk number 1: ex-qt-read-csv.Rnw:1-3
###################################################
## An example dialog to gather arguments for read.csv
require(qtbase)


###################################################
### code chunk number 2: ex-qt-read-csv.Rnw:22-36
###################################################
controls <- list()
controls$file <- Qt$QPushButton("click to select...")
##
controls$header <- Qt$QCheckBox()                 # no name
controls$header$setChecked(TRUE)
##
controls$sep <- Qt$QComboBox()
controls$sep$addItems(sprintf('%s', c(",", ";","","\t")))
controls$sep$setEditable(TRUE)
##               
controls$quote <- Qt$QLineEdit("\"'")
##
controls$fill <- Qt$QCheckBox()
controls$fill$setChecked(TRUE)


###################################################
### code chunk number 3: comment.char
###################################################
controls$comment.char <- Qt$QGroupBox()        # container
comment.char <- lapply(sprintf("%s", c("","#","%")),
                       Qt$QRadioButton, controls$comment.char)
comment.char[[1]]$setChecked(TRUE)
## manage
comment.char.bg <- Qt$QButtonGroup()
sapply(comment.char, comment.char.bg$addButton)
## layout
layout <- Qt$QVBoxLayout()
controls$comment.char$setLayout(layout)
sapply(comment.char, layout$addWidget)


###################################################
### code chunk number 4: ex-qt-read-csv.Rnw:66-70
###################################################
controls$name <- Qt$QLineEdit("")
controls$name$setPlaceholderText("Variable name to store data")
completer <- Qt$QCompleter(ls(.GlobalEnv))
controls$name$setCompleter(completer)


###################################################
### code chunk number 5: formLayout
###################################################
form_layout <- Qt$QFormLayout()
mapply(form_layout$addRow, names(controls), controls)


###################################################
### code chunk number 6: buttonBox
###################################################
button_box <- 
  Qt$QDialogButtonBox(Qt$QMessageBox$Cancel | 
                      Qt$QMessageBox$Ok)


###################################################
### code chunk number 7: windowLayout
###################################################
window <- Qt$QWidget()
window$windowTitle <- "Read csv file"
window$setLayout(window_layout <- Qt$QVBoxLayout())
window_layout$addLayout(form_layout)
window_layout$addWidget(button_box)


###################################################
### code chunk number 8: fileHandler
###################################################
filename <- NULL
qconnect(controls$file, "clicked", function() {
  name_filter <- "CSV file (*.csv);; All files (*.*)"
  filename <<- Qt$QFileDialog$getOpenFileName(window, 
               "Select a CSV file...", getwd(), name_filter)
  if(!is.null(filename))
    controls$file$setText(basename(filename))
})


###################################################
### code chunk number 9: buttonBox
###################################################
qconnect(button_box, "rejected", function() window$hide())
##
qconnect(button_box, "accepted", function() {
  if(!is.null(filename) && nchar(controls$name$text) > 0) {
    args <- list(file=filename,
                 header=controls$header$checked,
                 sep=controls$sep$currentText,
                 quote=controls$quote$text,
                 fill=controls$fill$checked
                 )
    args$comment.char <- comment.char.bg$checkedButton()$text
    ##
    val <- do.call("read.csv", args)
    assign(controls$name$text, val, .GlobalEnv)
    window$hide()
  } else {
    Qt$QMessageBox$warning(parent = window, 
       title = "Warning!",
       text = "You need to select a file and variable name")
  }
})
  


###################################################
### code chunk number 10: showAndRaise
###################################################
window$show()
window$raise()


