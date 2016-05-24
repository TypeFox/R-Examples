### R code from vignette source 'ch-Qt.Rnw'

###################################################
### code chunk number 1: ch-Qt.Rnw:4-34
###################################################
options(prompt=" ")
options(continue=" ")
options(width=80)
source("../booktabs.R")
##' find class that has method
##' @param obj object
##' @param meth_name if given prints out which class this comes from
##' @return a list with methods by class
findMeth <- function(obj, meth_name) {
  cls <- class(obj)
  ind <- sapply(cls, function(i) exists(i, envir=Qt))
  l <- lapply(cls[ind], function(i) ls(get(i, envir=Qt)()))
  names(l) <- cls[ind]
  for(i in 1:(length(l)-1))
    l[[i]] <- setdiff(l[[i]], l[[i+1]])

  if(!missing(meth_name)) {
    for(i in 1:length(l)) {
      if(meth_name %in% l[[i]])
        print(names(l)[i])
    }
  }
  invisible(l)
}

lookup <- function(obj, regxp) {
  x <- ls(obj)
  x[grepl(regxp, x)]
}
  


###################################################
### code chunk number 2: ex-qtbase.Rnw:6-7
###################################################
library(qtbase)


###################################################
### code chunk number 3: ex-qtbase.Rnw:25-29
###################################################
window <- Qt$QWidget()
label <- Qt$QLabel("Date:")
edit <- Qt$QLineEdit()
button <- Qt$QPushButton("Ok")


###################################################
### code chunk number 4: ex-qtbase.Rnw:48-49
###################################################
window$windowTitle <- "An example"


###################################################
### code chunk number 5: ex-qtbase.Rnw:66-67
###################################################
edit$setInputMask("0000-00-00")


###################################################
### code chunk number 6: ex-qtbase.Rnw:80-85
###################################################
layout <- Qt$QGridLayout()
layout$addWidget(label, row = 0, column = 0, 
                 rowSpan = 1, columnSpan = 1)
layout$addWidget(edit,   0, 1, 1, 1)
layout$addWidget(button, 1, 1, 1, 1)


###################################################
### code chunk number 7: ex-qtbase.Rnw:89-90
###################################################
window$setLayout(layout)


###################################################
### code chunk number 8: ex-qtbase.Rnw:95-96
###################################################
window$show()


###################################################
### code chunk number 9: ex-qtbase.Rnw:106-108
###################################################
handler <- function()  print(edit$text)
qconnect(button, "clicked", handler)


###################################################
### code chunk number 10: ex-qtbase.Rnw:132-136
###################################################
qsetClass("DateValidator", Qt$QValidator, 
          function(parent = NULL) {
            super(parent)
          })


###################################################
### code chunk number 11: ex-qtbase.Rnw:141-149
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
### code chunk number 12: ex-qtbase.Rnw:163-165
###################################################
validator <- DateValidator()
edit$setValidator(validator)


###################################################
### code chunk number 13: qtwidget-class
###################################################
Qt$QWidget


###################################################
### code chunk number 14: Overview.Rnw:104-105 (eval = FALSE)
###################################################
## head(names(Qt$QWidget), n = 3)


###################################################
### code chunk number 15: Overview.Rnw:112-113
###################################################
Qt$QWidget$DrawChildren


###################################################
### code chunk number 16: qwidget-object
###################################################
widget <- Qt$QWidget()


###################################################
### code chunk number 17: qwidget-class
###################################################
class(widget)


###################################################
### code chunk number 18: Overview.Rnw:137-138
###################################################
head(ls(widget), n=3)


###################################################
### code chunk number 19: windowTitle
###################################################
widget$windowTitle                           # initially NULL
widget$windowTitle <- "a new title"          # set property
widget$windowTitle


###################################################
### code chunk number 20: Overview.Rnw:161-162
###################################################
widget$show()


###################################################
### code chunk number 21: pushButton
###################################################
button <- Qt$QPushButton()


###################################################
### code chunk number 22: pushbutton-class
###################################################
is(button, "QWidget")
button$visible


###################################################
### code chunk number 23: Overview.Rnw:220-221
###################################################
button <- Qt$QPushButton()


###################################################
### code chunk number 24: constructor-parent
###################################################
widget <- Qt$QWidget()
button <- Qt$QPushButton(widget)


###################################################
### code chunk number 25: constructor-text
###################################################
button <- Qt$QPushButton("Button text")


###################################################
### code chunk number 26: constructor-icon-parent
###################################################
style <- Qt$QApplication$style()
icon <- style$standardIcon(Qt$QStyle$SP_DialogOkButton)
button <- Qt$QPushButton(icon, "Ok")


###################################################
### code chunk number 27: Overview.Rnw:263-266
###################################################
method_info <- qmethods(Qt$QPushButton)
dim(method_info)
head(method_info[,1:3], n = 3)


###################################################
### code chunk number 28: list-properties
###################################################
head(qproperties(button))


###################################################
### code chunk number 29: modal-property
###################################################
button$objectName <- "My button"
button$objectName
button$modal
cat(try(button$modal <- TRUE))


###################################################
### code chunk number 30: windowTitleMethod
###################################################
button$setObjectName("My button")


###################################################
### code chunk number 31: Overview.Rnw:322-325
###################################################
button <- Qt$QPushButton("click me")
qconnect(button, "clicked", function() message("ouch"))
button$show()


###################################################
### code chunk number 32: qt-signals-list
###################################################
tail(qsignals(Qt$QPushButton), n = 5)


###################################################
### code chunk number 33: qt-signals-qconnect-sig
###################################################
qconnect(button, "clicked", function(checked) print(checked))
qconnect(button, "clicked(bool)", 
         function(checked) print(checked))


###################################################
### code chunk number 34: qt-signals-disconnect
###################################################
proxy <- qconnect(button, "clicked", 
                  function() message("ouch"))
button$disconnect(proxy)


###################################################
### code chunk number 35: Enumeration-example
###################################################
Qt$Qt$AlignRight
Qt$QSizePolicy$Expanding


###################################################
### code chunk number 36: alignexample
###################################################
label <- Qt$QLabel("Our text")
label$alignment <- Qt$Qt$AlignRight | Qt$Qt$AlignTop


###################################################
### code chunk number 37: Overview.Rnw:428-429
###################################################
as.logical(label$alignment & Qt$Qt$AlignRight)


###################################################
### code chunk number 38: Overview.Rnw:459-460
###################################################
qsetClass("SubClass", Qt$QWidget)


###################################################
### code chunk number 39: qt-classes-show-class
###################################################
SubClass


###################################################
### code chunk number 40: qt-classes-constructor
###################################################
instance <- SubClass()


###################################################
### code chunk number 41: Overview.Rnw:483-488
###################################################
qsetClass("SubClass2", Qt$QWidget, 
          function(property, parent = NULL) {
            super(parent)
            this$property <- property
          })


###################################################
### code chunk number 42: Overview.Rnw:508-512
###################################################
qsetMethod("field", SubClass, function() field)
qsetMethod("setField", SubClass, function(value) {
  this$field <- value
})


###################################################
### code chunk number 43: Overview.Rnw:527-531
###################################################
qsetMethod("setVisible", SubClass, function(value) {
  message("Visible: ", value)
  super("setVisible", value)
})


###################################################
### code chunk number 44: subclass-show
###################################################
instance$show()


###################################################
### code chunk number 45: qt-overview-qsetSignal
###################################################
qsetSignal("somethingHappened", SubClass)


###################################################
### code chunk number 46: qt-overview-qsetSignal-sig
###################################################
qsetSignal("somethingHappenedAtIndex(int)", SubClass)


###################################################
### code chunk number 47: qt-overview-qsetSlot
###################################################
qsetSlot("doSomethingToIndex(int)", SubClass,function(index) {
  # ....
})


###################################################
### code chunk number 48: qt-overview-define-propE
###################################################
qsetProperty("property", SubClass)


###################################################
### code chunk number 49: qt-overview-define-prop-get
###################################################
instance <- SubClass()
instance$property                       # initially NULL
instance$property <- "value"
instance$property


###################################################
### code chunk number 50: qt-overview-define-prop-getter-setter
###################################################
qsetProperty("checkedProperty", SubClass, write=function(x) {
  if (!is(x, "character"))
    stop("'checkedProperty' must be a character vector")
  this$.checkedProperty <- x
})


###################################################
### code chunk number 51: qt-overview-define-prop-notify
###################################################
qsetSignal("propertyChanged", SubClass)
qsetProperty("property", SubClass, notify = "propertyChanged")


###################################################
### code chunk number 52: qt-overview-define-prop-typed
###################################################
qsetProperty("typedProperty", SubClass, type = "QString")


###################################################
### code chunk number 53: qt-overview-define-prop-typed-list
###################################################
tail(qproperties(SubClass()), 1)


###################################################
### code chunk number 54: ws_model
###################################################
library(qtbase)


###################################################
### code chunk number 55: ex-qt-ws-model.Rnw:16-20
###################################################
qsetClass("WSWatcher", Qt$QObject, function(parent = NULL) {
  super(parent)
  updateVariables()
})


###################################################
### code chunk number 56: ex-qt-ws-model.Rnw:26-27
###################################################
library(digest)


###################################################
### code chunk number 57: ex-qt-ws-model.Rnw:32-33
###################################################
qsetProperty("digests", WSWatcher)


###################################################
### code chunk number 58: objects_changed_property
###################################################
qsetSignal("objectsChanged", WSWatcher)


###################################################
### code chunk number 59: ex-qt-ws-model.Rnw:49-50
###################################################
qsetProperty("objects", WSWatcher, notify = "objectsChanged")


###################################################
### code chunk number 60: ex-qt-ws-model.Rnw:54-56
###################################################
qsetProperty("old_digests", WSWatcher)
qsetProperty("old_objects", WSWatcher)


###################################################
### code chunk number 61: update_variables
###################################################
qsetMethod("updateVariables", WSWatcher, function() {
  x <- sort(ls(envir = .GlobalEnv))
  objs <- sapply(mget(x, .GlobalEnv), digest)

  if((length(objs) != length(digests)) ||
     length(digests) == 0 ||
     any(objs != digests)) {
    this$old_digests <- digests         # old
    this$old_objects <- objects
    this$digests <- objs                # update cache
    this$objects <- x                   # emits signal         
  }
  invisible()
})


###################################################
### code chunk number 62: change_add
###################################################
qsetMethod("changedVariables", WSWatcher, function() {
  changed <- setdiff(old_digests, digests)
  old_objects[old_digests %in% changed]
})
##
qsetMethod("addedVariables", WSWatcher, function() {
  added <- setdiff(digests, old_digests)
  objects[digests %in% added]
})


###################################################
### code chunk number 63: addTaskCallback
###################################################
watcher <- WSWatcher()                          # an instance
addTaskCallback(function(expr, value, ok, visible) {
  watcher$updateVariables()
  TRUE
})


###################################################
### code chunk number 64: ex-qt-ws-model.Rnw:111-115 (eval = FALSE)
###################################################
## timer <- Qt$QTimer()
## timer$setSingleShot(FALSE)              # or TRUE for run once
## qconnect(timer,"timeout",function() watcher$updateVariables())
## timer$start(as.integer(3*1000))         # 3 seconds


###################################################
### code chunk number 65: connectSignal
###################################################
qconnect(watcher, "objectsChanged", function() 
         message("workspace objects were updated"))
new_object <- "The change should be announced"


###################################################
### code chunk number 66: qt-overview-visible
###################################################
widget <- Qt$QWidget()
widget$visible


###################################################
### code chunk number 67: qt-overview-show-hide
###################################################
widget$show()


###################################################
### code chunk number 68: Overview.Rnw:690-691
###################################################
widget$visible


###################################################
### code chunk number 69: Overview.Rnw:693-694
###################################################
widget$hide()


###################################################
### code chunk number 70: Overview.Rnw:696-697
###################################################
widget$visible


###################################################
### code chunk number 71: qt-overview-enabled
###################################################
button <- Qt$QPushButton("button")
button$enabled <- FALSE


###################################################
### code chunk number 72: qt-overview-size
###################################################
widget$size <- qsize(400, 400)
## or
widget$resize(400, 400)
widget$show()


###################################################
### code chunk number 73: Overview.Rnw:765-767
###################################################
font <- qfont(family = "helvetica", pointsize = 12, 
              weight = Qt$QFont$Bold, italic = TRUE)


###################################################
### code chunk number 74: Overview.Rnw:772-774
###################################################
label <- Qt$QLabel("Text for the label")
label$font <- font


###################################################
### code chunk number 75: qt-overview-stylesheets
###################################################
button <- Qt$QPushButton("Style sheet example")
button$show()
button$styleSheet <- 
  "QPushButton {color: red; background: white}"


###################################################
### code chunk number 76: Overview.Rnw:829-843 (eval = FALSE)
###################################################
## w <- Qt$QWidget()
## w$windowTitle <- "Using Style Sheets"
## lyt <- Qt$QHBoxLayout()
## w$setLayout(lyt)
## 
## b <- Qt$QPushButton("Style sheet example")
## lyt$addWidget(b)
## b1 <- Qt$QPushButton("Style sheet example")
## qsetStyleSheet(color = "red", background = "white", 
##                what = "QPushButton", widget = b1)
## lyt$addWidget(b1)
## 
## w$show()
## w$raise()


###################################################
### code chunk number 77: qt-overview-qsetStyleSheet (eval = FALSE)
###################################################
## qsetStyleSheet(color = "red", background = "white", 
##                what = "QPushButton", widget = button)


###################################################
### code chunk number 78: qt-overview-qsetStyleSheet-defaults
###################################################
qsetStyleSheet(color = "red", background = "white")


###################################################
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


###################################################
### code chunk number 84: has_loader
###################################################
has_uiloader <- FALSE#!is.null(Qt$QUiLoader)


###################################################
### code chunk number 85: qt-overview-designer-load (eval = FALSE)
###################################################
## loader <- Qt$QUiLoader()
## widget <- loader$load(Qt$QFile("textfinder.ui"))


###################################################
### code chunk number 86: Overview.Rnw:962-965
###################################################
if(has_uiloader) {
loader <- Qt$QUiLoader()
widget <- loader$load(Qt$QFile("textfinder.ui"))
}


###################################################
### code chunk number 87: qt-overview-designer-find (eval = FALSE)
###################################################
## find_button <- qfindChild(widget, "findButton") # by name
## line_edit <- qfindChild(widget, "lineEdit")


###################################################
### code chunk number 88: Overview.Rnw:986-989
###################################################
if(has_uiloader) {
find_button <- qfindChild(widget, "findButton") # by name
line_edit <- qfindChild(widget, "lineEdit")
}


###################################################
### code chunk number 89: qt-overview-designer-connect (eval = FALSE)
###################################################
## qconnect(find_button, "clicked", function() {
##   findText(line_edit$text)
## })


###################################################
### code chunk number 90: Overview.Rnw:999-1002
###################################################
if(has_uiloader) {
qconnect(find_button, "clicked", function() {
  findText(line_edit$text)
})
}


###################################################
### code chunk number 91: qt-overview-designer-custom-parent (eval = FALSE)
###################################################
## qsetClass("MyMainWindow", Qt$QWidget, function() {
##   loader <- Qt$QUiLoader()
##   widget <- loader$load(Qt$QFile("textfinder.ui"), this)
##   Qt$QMetaObject$connectSlotsByName(this)
## })


###################################################
### code chunk number 92: Overview.Rnw:1016-1019
###################################################
if(has_uiloader) {
qsetClass("MyMainWindow", Qt$QWidget, function() {
  loader <- Qt$QUiLoader()
  widget <- loader$load(Qt$QFile("textfinder.ui"), this)
  Qt$QMetaObject$connectSlotsByName(this)
})
}


###################################################
### code chunk number 93: qt-overview-designed-slot-handler (eval = FALSE)
###################################################
## qsetSlot("on_findButton_clicked", MyMainWindow, function() {
##   findText(line_edit$text)
## })


###################################################
### code chunk number 94: Overview.Rnw:1036-1039
###################################################
if(has_uiloader) {
qsetSlot("on_findButton_clicked", MyMainWindow, function() {
  findText(line_edit$text)
})
}


###################################################
### code chunk number 95: qt-overview-designed-custom-construct (eval = FALSE)
###################################################
## MyMainWindow()


###################################################
### code chunk number 96: Overview.Rnw:1048-1051
###################################################
if(has_uiloader) {
MyMainWindow()
}


###################################################
### code chunk number 97: ex-layout.Rnw:1-3
###################################################
## layout example
library(qtbase)


###################################################
### code chunk number 98: ex-layout.Rnw:22-26
###################################################
window <- Qt$QWidget()
window$setWindowTitle("Layout example")
grid_layout <- Qt$QGridLayout()
window$setLayout(grid_layout)


###################################################
### code chunk number 99: mainComponents
###################################################
fake_table <- Qt$QLabel("Table widget") 
notebook <- Qt$QTabWidget()
button_layout <- Qt$QHBoxLayout()


###################################################
### code chunk number 100: qt-layout-ex-add
###################################################
grid_layout$addWidget(fake_table, row=0, column=0, 
                     rowspan=1, colspan=1)
grid_layout$addWidget(notebook, 0, 1)
grid_layout$addLayout(button_layout, 1, 1)


###################################################
### code chunk number 101: layoutButtons
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
### code chunk number 102: nbLayout
###################################################
notebook_page <- Qt$QWidget()
notebook$addTab(notebook_page, "Tab label")
notebook$setTabToolTip(0, "A notebook page with a form")


###################################################
### code chunk number 103: nbFormLayout
###################################################
form_layout <- Qt$QFormLayout()
notebook_page$setLayout(form_layout)
l <- sapply(c("name", "rank", "snumber"),  Qt$QLineEdit)
form_layout$addRow("Name", l$name)
form_layout$addRow("Rank", l$rank)
form_layout$addRow("Serial number", l$snumber)


###################################################
### code chunk number 104: ex-layout.Rnw:96-98
###################################################
window$setMinimumSize(width=500, height=400)
window$show()


###################################################
### code chunk number 105: qt-layout-basic-constructor
###################################################
layout <- Qt$QHBoxLayout()


###################################################
### code chunk number 106: qt-layout-basic-setLayout
###################################################
widget <- Qt$QWidget()
widget$setLayout(layout)


###################################################
### code chunk number 107: qt-layout-basic-add
###################################################
layout$addWidget(Qt$QPushButton("Push Me"))


###################################################
### code chunk number 108: qt-layout-basic-item-at
###################################################
item <- layout$itemAt(0)


###################################################
### code chunk number 109: qt-layout-basic-widget
###################################################
button <- item$widget()


###################################################
### code chunk number 110: qt-layout-basic-remove
###################################################
layout$removeWidget(button)


###################################################
### code chunk number 111: qt-layout-basic-setParent
###################################################
button$setParent(NULL)


###################################################
### code chunk number 112: Layouts.Rnw:130-144
###################################################
DF <- rbind(
            c("Fixed","Require the size hint exactly"),
            c("Minimum", "Treat the size hint as the minimum, allowing expansion"),
            c("Maximum", "Treat the size hint as the maximum, allowing shrinkage"),
            c("Preferred", "Request the size hint, but allow for either expansion or shrinkage"),
            c("Expanding", "Treat like \\code{Preferred}, except the widget desires  as much space as possible"),
            c("MinimumExpanding", "Treat like \\code{Minimum}, except the widget desires  as much space as possible"),
            c("Ignored", "Ignore the size hint and request as much space as possible")
           )
colnames(DF) <- c("Policy","Meaning")
            cat(booktabs(DF,
                         colTypes=c("l","p{0.7\\textwidth}"),
                         caption="Possible size policies from \\class{QSizePolicy}",
                         label="tab:qt:size-policies"))


###################################################
### code chunk number 113: sizePolicy
###################################################
button <- Qt$QPushButton("No expansion")
button$setSizePolicy(vertical = Qt$QSizePolicy$Fixed, 
                     horizontal = Qt$QSizePolicy$Fixed)


###################################################
### code chunk number 114: Layouts.Rnw:173-174
###################################################
layout <- Qt$QHBoxLayout()


###################################################
### code chunk number 115: qt-layout-basic-alignment
###################################################
label <- Qt$QLabel("Label")
layout$addWidget(label)
layout$setAlignment(label, Qt$Qt$AlignLeft)


###################################################
### code chunk number 116: qt-layout-basic-spacing
###################################################
layout$spacing <- 5L


###################################################
### code chunk number 117: qt-layout-box
###################################################
hbox <- Qt$QHBoxLayout()
widget <- Qt$QWidget()
w$setLayout(hbox)


###################################################
### code chunk number 118: qt-layout-box-add
###################################################
buttons <- sapply(letters[1:3], Qt$QPushButton)
sapply(buttons, hbox$addWidget)


###################################################
### code chunk number 119: qt-layout-box-stretch
###################################################
hbox$setStretchFactor(buttons[[1]], 2.0) 


###################################################
### code chunk number 120: stretch_expand
###################################################
## Not shown
## Example of various combinations of stretch, alignment, sizePolicy
## no stretch, no expand
w <- Qt$QWidget(); w$setWindowTitle("No stretch, no expand")
w$setLayout(g <- Qt$QHBoxLayout())
buttons <- sapply(letters[1:3], Qt$QPushButton)
sapply(buttons, g$addWidget)
w$setMinimumSize(400, 50)
w$show()

## no stretch, expand for first
w <- Qt$QWidget();  w$setWindowTitle("No stretch, no expand, size Policy")
w$setLayout(g <- Qt$QHBoxLayout())

buttons <- sapply(letters[1:3], Qt$QPushButton)
sapply(buttons, g$addWidget)

b <- g$itemAt(0L)$widget()
b$setSizePolicy(Qt$QSizePolicy$Expanding, Qt$QSizePolicy$Fixed) 
for(i in 1:2) {
  b <- g$itemAt(i)$widget()
  b$setSizePolicy(Qt$QSizePolicy$Fixed, Qt$QSizePolicy$Fixed) 
}
w$setMinimumSize(400, 50)
w$show()                

## stretch
w <- Qt$QWidget(); w$setWindowTitle("Using stretch factors")
w$setLayout(g <- Qt$QHBoxLayout())
buttons <- sapply(letters[1:3], Qt$QPushButton)
for(i in 1:3) g$addWidget(buttons[[i]], stretch=i)


w$setMinimumSize(400, 50)
w$show()                

## no stretch; alignment
w <- Qt$QWidget(); w$setWindowTitle("Using alignment")
w$setLayout(g <- Qt$QHBoxLayout())
g$addWidget(Qt$QPushButton("NW"), stretch=0, Qt$Qt$AlignLeft | Qt$Qt$AlignTop)
g$addWidget(Qt$QPushButton("Center"), stretch=0, Qt$Qt$AlignHCenter | Qt$Qt$AlignVCenter)
g$addWidget(Qt$QPushButton("SE"), stretch=0, Qt$Qt$AlignRight | Qt$Qt$AlignBottom)
w$setMinimumSize(400, 400)
w$show()
##


###################################################
### code chunk number 121: Layouts.Rnw:291-299
###################################################
## Not shown
## example of inserting "A" at index=1. 0 1 2 3 -> 0 A 1 2 3
w <- Qt$QWidget()
w$setWindowTitle("Insert example")
w$setLayout(g <- Qt$QHBoxLayout())
for(i in 0:2) g$addWidget(Qt$QPushButton(as.character(i)))
g$insertWidget(1L, Qt$QPushButton("inserted"))
w$show()


###################################################
### code chunk number 122: create_hb
###################################################
hbox <- Qt$QHBoxLayout()


###################################################
### code chunk number 123: qt-layout-box-spacing
###################################################
hbox$addSpacing(12L)
hbox$addWidget(Qt$QPushButton("d"))


###################################################
### code chunk number 124: qt-layout-box-stretch
###################################################
hbox$addStretch(2)
hbox$addWidget(Qt$QPushButton("Help..."))


###################################################
### code chunk number 125: qt-layout-box-strut
###################################################
hbox$addStrut(10)                    # at least 10 pixels high


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


###################################################
### code chunk number 134: FormLayoutExample
###################################################
window <- Qt$QWidget()
window$setWindowTitle("Wrapper for 'dnorm' function")
window$setLayout(layout <- Qt$QFormLayout())
sapply(c("quantile", "mean", "sd"), function(statistic) {
  layout$addRow(statistic, Qt$QLineEdit())
})
layout$addRow(Qt$QCheckBox("log"))


###################################################
### code chunk number 135: Layouts.Rnw:492-493
###################################################
window$show(); window$raise()


###################################################
### code chunk number 136: qt-layout-form-at
###################################################
item <- layout$itemAt(0, Qt$QFormLayout$FieldRole)
quantile_edit <- item$widget()


###################################################
### code chunk number 137: groupBoxExample
###################################################
## Not shown, Group box example of title, alignment, flat, checkable
f <- Qt$QGroupBox("Example group box")
lyt <- Qt$QVBoxLayout(); f$setLayout(lyt)
lyt$addWidget(changeTitle <- Qt$QPushButton("Change title"))
lyt$addWidget(changeAlignment <- Qt$QPushButton("Cycle Alignment"))
lyt$addWidget(toggleFlat <- Qt$QPushButton("Toggle flat"))
lyt$addWidget(toggleCheckable <- Qt$QPushButton("Toggle checkable"))
f$show()

qconnect(changeTitle, "clicked", function(checked) {
  f$setTitle("New title")
})
qconnect(changeAlignment, "clicked", function(checked) {
  aligns <- c(Qt$Qt$AlignLeft, Qt$Qt$AlignHCenter, Qt$Qt$AlignRight)
  curAlign <- f$alignment
  ind <-   which(curAlign == aligns)
  f$setAlignment(aligns[c(2,3,1)[ind]])
})
qconnect(toggleFlat, "clicked", function(checked) {
  f$setFlat(!f$flat)
})
qconnect(toggleCheckable, "clicked", function(checked) {
  f$setCheckable(!f$checkable)
})


###################################################
### code chunk number 138: qt-widget-separator
###################################################
separator <- Qt$QFrame()
separator$frameShape <- Qt$QFrame$HLine


###################################################
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


###################################################
### code chunk number 145: HelpBrowser
###################################################
qsetClass("HelpBrowser", Qt$QTabWidget, function(parent=NULL){
  super(parent)
  #
  this$tabsClosable <- TRUE
  qconnect(this, "tabCloseRequested", function(index) {
    this$removeTab(index)
  })
  this$movable <- TRUE; this$usesScrollButtons <- TRUE
  #
  this$browser <- getOption("browser")
  options("browser" =  function(url) openPage(url))
})


###################################################
### code chunk number 146: openPage
###################################################
qsetMethod("openPage", HelpBrowser, function(url) {
  tokens <- strsplit(url, "/")[[1]]
  tab_title <- sprintf("%s: %s", tokens[length(tokens)-2], 
                       tokens[length(tokens)])
  webview <- Qt$QWebView()
  webview$setUrl(Qt$QUrl(url))
  this$currentIndex <- addTab(webview, tab_title)
})


###################################################
### code chunk number 147: illustrateHelpBrowser
###################################################
help_browser <- HelpBrowser()
help_browser$windowTitle <- "Help Browser example"
help_browser$show()
help_browser$raise()
##
options("help_type"="html")
help("mean")
help("boxplot")


###################################################
### code chunk number 148: Layouts.Rnw:707-711
###################################################
image <- Qt$QLabel()
image$pixmap <- Qt$QPixmap("someimage.png")
scroll_area <- Qt$QScrollArea()
scroll_area$setWidget(image)


###################################################
### code chunk number 149: qt-layout-scrolled-zoom
###################################################
zoomImage <- function(x = 2.0) {
  image$resize(x * image$pixmap$size())
  updateScrollBar <- function(sb) {
    sb$value <- x * sb$value + (x - 1) * sb$pageStep / 2
  }
  updateScrollBar(scroll_area$horizontalScrollBar())
  updateScrollBar(scroll_area$verticalScrollBar())
}


###################################################
### code chunk number 150: Layouts.Rnw:741-743
###################################################
scroll_area$show()
scroll_area$raise()


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


###################################################
### code chunk number 155: qt-dialogs-static-warning (eval = FALSE)
###################################################
## response <- Qt$QMessageBox$warning(parent = NULL, 
##               title = "Warning!", text = "Warning message...")


###################################################
### code chunk number 156: qt-dialogs-construct
###################################################
dialog <- Qt$QMessageBox(icon = Qt$QMessageBox$Warning,
                         title = "Warning!",
                         text = "Warning text...",
                         buttons = Qt$QMessageBox$Ok,
                         parent = NULL)


###################################################
### code chunk number 157: qt-dialogs-extra-text
###################################################
dialog$informativeText <- "Less important warning information"
dialog$detailedText <- "Extra details most do not care to see"


###################################################
### code chunk number 158: qt-dialogs-exec
###################################################
dialog$exec()                       # returns response code


###################################################
### code chunk number 159: qt-dialogs-listen
###################################################
qconnect(dialog, "finished", function(response) {
  dialog$close()
})


###################################################
### code chunk number 160: qt-dialogs-show
###################################################
dialog$show()
dialog$raise()
dialog$activateWindow()


###################################################
### code chunk number 161: QMEssageBoxAPI
###################################################
dialog <- Qt$QMessageBox()
dialog$windowTitle <- "[This space for rent]"
dialog$text <- "This is the main text"
dialog$informativeText <- "This should give extra info"
dialog$detailedText <- "And this provides\neven more detail"

dialog$icon <- Qt$QMessageBox$Critical
dialog$standardButtons <- 
  Qt$QMessageBox$Cancel | Qt$QMessageBox$Ok
## 'Cancel' instead of 'Ok' is the default
dialog$setDefaultButton(Qt$QMessageBox$Cancel)
##
if(dialog$exec() == Qt$QMessageBox$Ok) 
  print("A Ok")


###################################################
### code chunk number 162: qt-dialogs-input-get-text
###################################################
text <- Qt$QInputDialog$getText(parent = NULL, 
                        title = "Gather text",
                        label = "Enter some text")


###################################################
### code chunk number 163: qt-dialogs-input-get-range
###################################################
even_integer <- Qt$QInputDialog$getInt(parent = NULL, 
                       title="Gather integer",
                       label="Enter an integer from 1 to 10",
                       value=0, min = 2, max = 10, step = 2)


###################################################
### code chunk number 164: qt-dialogs-input-get-item
###################################################
item <- Qt$QInputDialog$getItem(parent = NULL, 
                        title = "Select item",
                        label = "Select a letter",
                        items = LETTERS, current = 17)


###################################################
### code chunk number 165: qt-dialogs-input-explicit
###################################################
dialog <- Qt$QInputDialog()
dialog$setWindowTitle("Select item")
dialog$setLabelText("Select a letter")
dialog$setComboBoxItems(LETTERS)
dialog$setTextValue(LETTERS[18])
dialog$setOptions(Qt$QInputDialog$UseListViewForComboBoxItems)


###################################################
### code chunk number 166: Dialogs.Rnw:245-247
###################################################
if (dialog$exec())
  print(dialog$textValue())


###################################################
### code chunk number 167: qt-widgets-buttonbox
###################################################
btn_box <- Qt$QDialogButtonBox(Qt$QDialogButtonBox$Ok | 
                               Qt$QDialogButtonBox$Cancel | 
                               Qt$QDialogButtonBox$Help)


###################################################
### code chunk number 168: qt-widgets-buttonbox-signals
###################################################
qconnect(btn_box, "accepted", function() message("accepted"))
qconnect(btn_box, "rejected", function() message("rejected"))
qconnect(btn_box, "helpRequested", function() message("help"))
qconnect(btn_box, "clicked", 
         function(button) message(button$text))


###################################################
### code chunk number 169: Dialogs.Rnw:308-309
###################################################
btn_box$show(); btn_box$raise()


###################################################
### code chunk number 170: extend-date-dialog
###################################################
qsetClass("DateDialog", Qt$QDialog, 
          function(parent = NULL) {
            super(parent=parent)
            setWindowTitle("Choose a date")
            this$calendar <- Qt$QCalendarWidget()
            #
            btn_box <- 
              Qt$QDialogButtonBox(Qt$QMessageBox$Cancel | 
                                  Qt$QMessageBox$Ok)
            qconnect(btn_box, "accepted", function() {
              this$close()
              this$setResult(Qt$QMessageBox$Ok)    
            })
            qconnect(btn_box, "rejected", 
                     function() this$close())
            #
            layout <- Qt$QVBoxLayout()
            sapply(list(calendar, btn_box), layout$addWidget)
            setLayout(layout)
          })


###################################################
### code chunk number 171: extend-date-get
###################################################
qsetMethod("selectedDate", DateDialog, 
           function(x) calendar$selectedDate$toString())


###################################################
### code chunk number 172: extend-date-exec
###################################################
date_dialog <- DateDialog()
if (date_dialog$exec())
  message(date_dialog$selectedDate())


###################################################
### code chunk number 173: Wizard
###################################################
wizard <- Qt$QWizard()
wizard$setWindowTitle("A wizard")


###################################################
### code chunk number 174: Dialogs.Rnw:401-407
###################################################
get_age_page <- Qt$QWizardPage(wizard)
get_age_page$setTitle("What is your age?")
layout <- Qt$QFormLayout()
get_age_page$setLayout(layout)
layout$addRow("Age", (age <- Qt$QLineEdit()))
wizard$addPage(get_age_page)


###################################################
### code chunk number 175: Dialogs.Rnw:411-424
###################################################
get_toys_page <- Qt$QWizardPage(wizard)
get_toys_page$setTitle("What toys do you like?")
layout <- Qt$QFormLayout()
get_toys_page$setLayout(layout)
layout$addRow("Toys", (toys <- Qt$QLineEdit()))
wizard$addPage(get_toys_page)
##
get_games_page <- Qt$QWizardPage(wizard)
get_games_page$setTitle("What games do you like?")
layout <- Qt$QFormLayout()
get_games_page$setLayout(layout)
layout$addRow("Games", (games <- Qt$QLineEdit()))
wizard$addPage(get_games_page)


###################################################
### code chunk number 176: Dialogs.Rnw:429-432
###################################################
response <- wizard$exec()
if(response)
  message(toys$text)


###################################################
### code chunk number 177: QFIleDialog (eval = FALSE)
###################################################
## filename <- Qt$QFileDialog$getOpenFileName(NULL, 
##                           "Open a file...", getwd())


###################################################
### code chunk number 178: Dialogs.Rnw:465-467 (eval = FALSE)
###################################################
## filenames <- Qt$QFileDialog$getOpenFileNames(NULL, 
##                         "Open file(s)...", getwd())


###################################################
### code chunk number 179: Dialogs.Rnw:471-473 (eval = FALSE)
###################################################
## filename <- Qt$QFileDialog$getSaveFileName(NULL,
##                         "Save as...", getwd())


###################################################
### code chunk number 180: Dialogs.Rnw:477-479 (eval = FALSE)
###################################################
## dirname <- Qt$QFileDialog$getExistingDirectory(NULL,
##                         "Select directory", getwd())


###################################################
### code chunk number 181: Dialogs.Rnw:487-494
###################################################
name_filter <- paste("R files (*.R .RData)",
                    "Sweave files (*.Rnw)",
                    "All files (*.*)", 
                    sep=";;")
##
filenames <- Qt$QFileDialog$getOpenFileNames(NULL, 
             "Open file(s)...", getwd(), name_filter)


###################################################
### code chunk number 182: Dialogs.Rnw:507-508
###################################################
dirname <- ""


###################################################
### code chunk number 183: QFileDialogAPI
###################################################
dialog <- Qt$QFileDialog(NULL, "Choose an R file", getwd(), 
                      name_filter)
dialog$fileMode <- Qt$QFileDialog$ExistingFiles
dialog$setHistory(dirname)


###################################################
### code chunk number 184: Dialogs.Rnw:520-522 (eval = FALSE)
###################################################
## if(dialog$exec())
##   print(dialog$selectedFiles())


###################################################
### code chunk number 185: qt-widget-label
###################################################
label <- Qt$QLabel("<font color='red'>Red</font>")


###################################################
### code chunk number 186: Widgets.Rnw:45-47
###################################################
## labels can use rich text
label$show(); label$raise()


###################################################
### code chunk number 187: qt-widget-button
###################################################
button <- Qt$QPushButton("Ok")


###################################################
### code chunk number 188: Widgets.Rnw:71-72
###################################################
button$show(); button$raise()


###################################################
### code chunk number 189: qt-widget-button-disable
###################################################
button$enabled <- FALSE


###################################################
### code chunk number 190: qt-widget-button-connect
###################################################
qconnect(button, "clicked", function() message("Ok clicked") )


###################################################
### code chunk number 191: qt-widget-button-icon
###################################################
icon_file <- system.file("images/ok.gif", package="gWidgets")
button$icon <- Qt$QIcon(icon_file)


###################################################
### code chunk number 192: qt-widget-button-icon-from-style
###################################################
style <- Qt$QApplication$style()
button$icon <- style$standardIcon(Qt$QStyle$SP_DialogOkButton)


###################################################
### code chunk number 193: Widgets.Rnw:138-140
###################################################
require(grid)
require(ggplot2)


###################################################
### code chunk number 194: Widgets.Rnw:142-146 (eval = FALSE)
###################################################
## require(qtutils)
## device <- QT()
## grid:::grid.newpage()
## grid:::grid.draw(ggplot2:::GeomHistogram$icon())


###################################################
### code chunk number 195: Widgets.Rnw:152-158 (eval = FALSE)
###################################################
## pixmap <- Qt$QPixmap(device$size$toSize())
## pixmap$fill()
## painter <- Qt$QPainter()
## painter$begin(pixmap)
## device$render(painter)
## painter$end()


###################################################
### code chunk number 196: Widgets.Rnw:163-165 (eval = FALSE)
###################################################
## button <- Qt$QPushButton("Histogram")
## button$setIcon(Qt$QIcon(pixmap))


###################################################
### code chunk number 197: Widgets.Rnw:168-170 (eval = FALSE)
###################################################
## button$show()
## button$raise()


###################################################
### code chunk number 198: qt-widget-checkbox
###################################################
checkbox <- Qt$QCheckBox("Option")


###################################################
### code chunk number 199: qt-widget-checkbox-checked
###################################################
checkbox$checked


###################################################
### code chunk number 200: qt-widget-checkbox-state-changed
###################################################
qconnect(checkbox, "stateChanged", function(state) {
  if (state == Qt$Qt$Checked)
    message("checked")
})


###################################################
### code chunk number 201: Widgets.Rnw:240-244
###################################################
window <- Qt$QWidget()
group_box <- Qt$QGroupBox("Cylinders:")
layout <- Qt$QVBoxLayout()
window$setLayout(layout)


###################################################
### code chunk number 202: Widgets.Rnw:248-250
###################################################
btn_group <- Qt$QButtonGroup()
btn_group$exclusive <- FALSE


###################################################
### code chunk number 203: Widgets.Rnw:259-268
###################################################
data(Cars93, package="MASS")
cylinders <- levels(Cars93$Cylinders)
sapply(seq_along(cylinders), function(i) {
  button <- Qt$QCheckBox(sprintf("%s Cylinders", cylinders[i]))
  layout$addWidget(button)
  btn_group$addButton(button, i)
})
sapply(btn_group$buttons(), 
       function(button) button$checked <- TRUE)


###################################################
### code chunk number 204: Widgets.Rnw:280-286
###################################################
checked <- sapply(btn_group$buttons(), function(i) i$checked)
if(any(checked)) {
  checked_cyls <- Cars93$Cylinders %in% cylinders[checked]
  message(sprintf("You've selected %d cases", 
                  sum(checked_cyls)))
}


###################################################
### code chunk number 205: Widgets.Rnw:299-305
###################################################
qconnect(btn_group, "buttonClicked(QAbstractButton*)", 
         function(button) {
           msg <- sprintf("Level '%s': %s", 
                          button$text, button$checked)
           message(msg)
})


###################################################
### code chunk number 206: Widgets.Rnw:307-309
###################################################
window$show()
window$raise()


###################################################
### code chunk number 207: RadioWithList
###################################################
window <- Qt$QGroupBox("Weight:")
radio_buttons <- 
  list(Qt$QRadioButton("Weight < 3000", w),
       Qt$QRadioButton("3000 <= Weight < 4000", w),
       Qt$QRadioButton("4000 <= Weight", w))


###################################################
### code chunk number 208: qt-widget-radio-layout
###################################################
layout <- Qt$QVBoxLayout()
window$setLayout(layout)
sapply(radio_buttons, layout$addWidget)
radio_buttons[[1]]$setChecked(TRUE)


###################################################
### code chunk number 209: qt-widget-radio-checked
###################################################
radio_buttons[[1]]$checked


###################################################
### code chunk number 210: Widgets.Rnw:348-355
###################################################
sapply(radio_buttons, function(button) {
  qconnect(button, "toggled", function(checked) {
    if(checked) {
      message(sprintf("You checked %s.", button$text))
    }
  })
})


###################################################
### code chunk number 211: qt-widget-radio-group
###################################################
btn_group <- Qt$QButtonGroup()
lapply(radio_buttons, btn_group$addButton)


###################################################
### code chunk number 212: qt-widget-radio-group-checked
###################################################
btn_group$checkedButton()$text


###################################################
### code chunk number 213: Widgets.Rnw:379-381
###################################################
window$show()
window$raise()


###################################################
### code chunk number 214: makeFigureFromPrevious
###################################################
## Code used to make figure, just combines previous
make_cb <- function() {
  w <- Qt$QGroupBox("Cylinders:")
  layout <- Qt$QVBoxLayout()
  w$setLayout(layout)
  btn_group <- Qt$QButtonGroup()
  btn_group$exclusive <- FALSE
  data(Cars93, package="MASS")
  cylinders <- levels(Cars93$Cylinders)
  sapply(seq_along(cylinders), function(i) {
    button <- Qt$QCheckBox(sprintf("%s Cylinders", cylinders[i]))
    layout$addWidget(button)
    btn_group$addButton(button, i)
  })
  sapply(btn_group$buttons(), function(i) i$checked <- TRUE)
  w
}

make_rb <- function() {
  w <- Qt$QGroupBox("Weight:")
  l <- list(Qt$QRadioButton("Weight < 3000", w),
            Qt$QRadioButton("3000 <= Weight < 4000", w),
            Qt$QRadioButton("4000 <= Weight", w))
  layout <- Qt$QVBoxLayout()
  w$setLayout(layout)
  sapply(l, function(i) layout$addWidget(i))
  l[[1]]$setChecked(TRUE)
  l[[1]]$checked
  sapply(l, function(i) {
    qconnect(i, "toggled", function(checked) {
      if(checked) {
        message(sprintf("You checked %s.", i$text))
      }
    })
  })
  btn_group <- Qt$QButtonGroup()
  lapply(l, btn_group$addButton)
  w
}

w <- Qt$QWidget()
layout <- Qt$QHBoxLayout()
layout$addWidget(make_cb())
layout$addWidget(make_rb())
w$setLayout(layout) 
w$windowTitle <- "Groups of checkboxes and radio buttons"


###################################################
### code chunk number 215: Widgets.Rnw:456-460
###################################################
df <- data.frame(name=state.name, region=state.region,
                 population=state.x77[,'Population'], 
                 stringsAsFactors=FALSE)
states_by_region <- split(df, df$region)


###################################################
### code chunk number 216: QComboBox
###################################################
state_combo <- Qt$QComboBox()
region_combo <- Qt$QComboBox()
region_combo$addItems(names(states_by_region))


###################################################
### code chunk number 217: qt-widget-combo-currentitem
###################################################
region_combo$currentText
region_combo$currentIndex                     # 0-based


###################################################
### code chunk number 218: qt-widget-combo-clear-index
###################################################
region_combo$currentIndex <- -1


###################################################
### code chunk number 219: Widgets.Rnw:491-495
###################################################
qconnect(region_combo, "activated(int)", function(index) {
  state_combo$clear()
  state_combo$addItems(states_by_region[[index+1]]$name)
})


###################################################
### code chunk number 220: Widgets.Rnw:502-509
###################################################
window <- Qt$QGroupBox("Two combo boxes")
layout <- Qt$QFormLayout()
window$setLayout(layout)
layout$addRow("Region:", region_combo)
layout$addRow("State:", state_combo)
layout$fieldGrowthPolicy <-  # grow combo boxes
  Qt$QFormLayout$AllNonFixedFieldsGrow


###################################################
### code chunk number 221: Widgets.Rnw:512-513 (eval = FALSE)
###################################################
## window$show(); window$raise()


###################################################
### code chunk number 222: qt-widget-slider
###################################################
slider <- Qt$QSlider()
slider$minimum <- 0
slider$maximum <- 100


###################################################
### code chunk number 223: qt-widget-slider-step
###################################################
slider$singleStep <- 1
slider$pageStep <- 5


###################################################
### code chunk number 224: qt-widget-slider-value
###################################################
slider$value
slider$value <- 50


###################################################
### code chunk number 225: qt-widget-slider-aesthetics
###################################################
slider$orientation <- Qt$Qt$Horizontal
slider$tickPosition <- Qt$QSlider$TicksBelow
slider$tickInterval <- 10


###################################################
### code chunk number 226: Widgets.Rnw:581-585
###################################################
spinbox <- Qt$QSpinBox()
spinbox$minimum <- slider$minimum
spinbox$maximum <- slider$maximum
spinbox$singleStep <- slider$singleStep


###################################################
### code chunk number 227: qt-widget-spin-suffix
###################################################
spinbox$suffix <- "%"


###################################################
### code chunk number 228: Widgets.Rnw:601-604
###################################################
f <- function(value, obj) obj$value <- value
qconnect(spinbox, "valueChanged", f, user.data = slider)
qconnect(slider, "valueChanged", f, user.data = spinbox)


###################################################
### code chunk number 229: SliderSpinButton
###################################################
w <- Qt$QWidget()
layout <- Qt$QHBoxLayout()
w$setLayout(layout)


###################################################
### code chunk number 230: Widgets.Rnw:618-624
###################################################
## not shown
layout$addWidget(slider)
layout$addWidget(spinbox)

w$show()
w$raise()


###################################################
### code chunk number 231: qt-widget-lineedit
###################################################
edit <- Qt$QLineEdit("Initial contents")


###################################################
### code chunk number 232: qt-widget-lineedit-text
###################################################
edit$text


###################################################
### code chunk number 233: qt-widget-lineedit-select
###################################################
edit$setSelection(start = 0, length = nchar(edit$text))


###################################################
### code chunk number 234: qt-widget-linedit-selectedText
###################################################
edit$selectedText


###################################################
### code chunk number 235: qt-widget-lineedit-placeholder
###################################################
edit$text <- ""
edit$setPlaceholderText("Enter some text here")


###################################################
### code chunk number 236: qt-widget-lineedit-editingFinished
###################################################
qconnect(edit, "editingFinished", function() {
  message("Entered text: '", edit$text, "'")
})


###################################################
### code chunk number 237: CompleterExample
###################################################
class_browser <- Qt$QWidget()
layout <- Qt$QFormLayout()
class_browser$setLayout(layout)

layout$addRow("Class name", class_edit <- Qt$QLineEdit())
layout$addRow("Method name", method_edit <- Qt$QLineEdit())


###################################################
### code chunk number 238: Widgets.Rnw:722-724
###################################################
class_completer <- Qt$QCompleter(ls(Qt))
class_edit$setCompleter(class_completer)


###################################################
### code chunk number 239: Widgets.Rnw:730-739
###################################################
qconnect(class_edit, "editingFinished", function() {
  class_name <- class_edit$text
  if(class_name == "") return()
  class_object <- get(class_name, envir = Qt)
  if(!is.null(class_object)) {
    method_completer <- Qt$QCompleter(ls(class_object()))
    method_edit$setCompleter(method_completer)
  }
})


###################################################
### code chunk number 240: Widgets.Rnw:743-748
###################################################
w <- Qt$QWidget()
w$windowTitle <- "Completion example"
w$setLayout(layout)
w$show()
w$raise()


###################################################
### code chunk number 241: qt-widget-lineedit-mask
###################################################
edit$inputMask <- "999-99-9999"


###################################################
### code chunk number 242: ex-qt-read-csv.Rnw:1-3
###################################################
## An example dialog to gather arguments for read.csv
require(qtbase)


###################################################
### code chunk number 243: ex-qt-read-csv.Rnw:22-36
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
### code chunk number 244: comment.char
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
### code chunk number 245: ex-qt-read-csv.Rnw:66-70
###################################################
controls$name <- Qt$QLineEdit("")
controls$name$setPlaceholderText("Variable name to store data")
completer <- Qt$QCompleter(ls(.GlobalEnv))
controls$name$setCompleter(completer)


###################################################
### code chunk number 246: formLayout
###################################################
form_layout <- Qt$QFormLayout()
mapply(form_layout$addRow, names(controls), controls)


###################################################
### code chunk number 247: buttonBox
###################################################
button_box <- 
  Qt$QDialogButtonBox(Qt$QMessageBox$Cancel | 
                      Qt$QMessageBox$Ok)


###################################################
### code chunk number 248: windowLayout
###################################################
window <- Qt$QWidget()
window$windowTitle <- "Read csv file"
window$setLayout(window_layout <- Qt$QVBoxLayout())
window_layout$addLayout(form_layout)
window_layout$addWidget(button_box)


###################################################
### code chunk number 249: fileHandler
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
### code chunk number 250: buttonBox
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
### code chunk number 251: showAndRaise
###################################################
window$show()
window$raise()


###################################################
### code chunk number 252: qt-widget-webview
###################################################
webview <- Qt$QWebView()
webview$load(Qt$QUrl("http://www.r-project.org"))


###################################################
### code chunk number 253: qt-widget-webview-plugin-html
###################################################
html <- readLines(out <- textConnection("
<html xmlns='http://www.w3.org/1999/xhtml'>
  <body>
    <h1>Qt class browser embedded into a QWebView</h1>
    Search for a class:<br/>
    <object type='application/x-qt-class-browser' width='500'
            height='100'/>
  </body>
</html>
")); close(out)
html <- paste(html, collapse = "\n")


###################################################
### code chunk number 254: qt-widget-webview-plugin-factory
###################################################
qsetClass("RPluginFactory", Qt$QWebPluginFactory)


###################################################
### code chunk number 255: qt-widget-webview-plugin-factory-plugins
###################################################
qsetMethod("plugins", RPluginFactory, function() {
  plugin <- Qt$QWebPluginFactory$Plugin()
  plugin$setName("Class Browser")
  mimeType <- Qt$QWebPluginFactory$MimeType()
  mimeType$setName("application/x-qt-class-browser")
  plugin$setMimeTypes(list(mimeType))
  list(plugin)
})


###################################################
### code chunk number 256: qt-widget-webview-plugin-factory-create
###################################################
qsetMethod("create", RPluginFactory, 
           function(mime_type, url, arg_names, arg_vals) {
             if (mime_type== "application/x-qt-class-browser")
              class_browser
             else Qt$QWidget()
           })


###################################################
### code chunk number 257: qt-widget-webview-plugin-factory-register
###################################################
globalSettings <- Qt$QWebSettings$globalSettings()
globalSettings$setAttribute(Qt$QWebSettings$PluginsEnabled, 
                            TRUE)
webview$page()$setPluginFactory(RPluginFactory())
webview$setHtml(html)


###################################################
### code chunk number 258: qt-widget-graphics-device (eval = FALSE)
###################################################
## library(qtutils)
## qt_device <- QT()
## plot(mpg ~ hp, data = mtcars)


###################################################
### code chunk number 259: qt-widget-graphics-device-notebook (eval = FALSE)
###################################################
## notebook <- Qt$QTabWidget()
## notebook$addTab(qt_device, "Plot 1")
## print(notebook)


###################################################
### code chunk number 260: qt-widget-graphics-device-opengl (eval = FALSE)
###################################################
## qt_opengl_device <- QT(opengl = TRUE)


###################################################
### code chunk number 261: DragConstructor
###################################################
qsetClass("DragLabel", Qt$QLabel, 
          function(text = "", parent = NULL) {
            super(parent)
            setText(text)
            ##
            setAlignment(Qt$Qt$AlignCenter)
            setMinimumSize(200, 200)
          })


###################################################
### code chunk number 262: drag-mouse-press-event
###################################################
qsetMethod("mousePressEvent", DragLabel, function(event) {
  mime_data <- Qt$QMimeData()
  mime_data$setText(text)

  drag <- Qt$QDrag(this)
  drag$setMimeData(mime_data)

  drag$exec()
})


###################################################
### code chunk number 263: DropConstructor
###################################################
qsetClass("DropLabel", Qt$QLabel, 
          function(text="", parent=NULL) {
            super(parent)
            
            setText(text)
            this$acceptDrops <- TRUE
            
            this$bgrole <- backgroundRole()
            this$alignment <- Qt$Qt$AlignCenter
            setMinimumSize(200, 200)
            this$autoFillBackground <- TRUE
            clear()
          })


###################################################
### code chunk number 264: qt-dnd-utility
###################################################
qsetMethod("clear", DropLabel, function() {
  setText(this$orig_text)
  setBackgroundRole(this$bgrole)
})
qsetMethod("setText", DropLabel, function(text) {
  this$orig_text <- text
  super("setText", text)                 # next method
})


###################################################
### code chunk number 265: Widgets.Rnw:1033-1043
###################################################
qsetMethod("dragEnterEvent", DropLabel, function(event) {
  mime_data <- event$mimeData()
  if (mime_data$hasImage() || mime_data$hasHtml() | 
      mime_data$hasText()) 
  {
    super("setText", "<Drop Text Here>")
    setBackgroundRole(Qt$QPalette$Highlight)
    event$acceptProposedAction()
  }
})


###################################################
### code chunk number 266: Widgets.Rnw:1051-1054
###################################################
qsetMethod("dragLeaveEvent", DropLabel, function(event) {
  clear()
})


###################################################
### code chunk number 267: dropevent
###################################################
qsetMethod("dropEvent", DropLabel, function(event) {
  mime_data <- event$mimeData()
  
  if(mime_data$hasImage()) {
    setPixmap(mime_data$imageData())
  } else if(mime_data$hasHtml()) {
    setText(mime_data$html)
    setTextFormat(Qt$Qt$RichText)
  } else if(mime_data$hasText()) {
    setText(mime_data$text())
    setTextFormat(Qt$Qt$PlainText)
  } else {
    setText("No match")                 # replace...
  }

  setBackgroundRole(this$bgrole)
  event$acceptProposedAction()
})


###################################################
### code chunk number 268: Widgets-MVC.Rnw:2-4
###################################################
require(qtbase)
require(MASS)


###################################################
### code chunk number 269: qt-mvc-qdfm
###################################################
model <- qdataFrameModel(mtcars)
view <- Qt$QTableView()
view$setModel(model)


###################################################
### code chunk number 270: qt-mvc-qdfm-access
###################################################
DF <- qdataFrame(model)
DF[1:3, 1:10]


###################################################
### code chunk number 271: Widgets-MVC.Rnw:72-73
###################################################
qdataFrame(model)$hpToMpg <- with(qdataFrame(model), hp / mpg)


###################################################
### code chunk number 272: qt-mvc-table-header-align
###################################################
header <- view$horizontalHeader()
header$defaultAlignment <- Qt$Qt$AlignLeft


###################################################
### code chunk number 273: qt-mvc-table-alternating
###################################################
view$alternatingRowColors <- TRUE


###################################################
### code chunk number 274: qt-mvc-table-gc
###################################################
broken_view <- Qt$QTableView()
broken_view$setModel(qdataFrameModel(mtcars))
gc()


###################################################
### code chunk number 275: qt-mvc-table-gc-result
###################################################
broken_view$model()              # NULL, garbage collected


###################################################
### code chunk number 276: qt-mvc-table-parent
###################################################
parental_view <- Qt$QTableView()
broken_view$setModel(qdataFrameModel(mtcars, 
                                     parent = parental_view))
gc()


###################################################
### code chunk number 277: qt-mvc-table-parent-result
###################################################
broken_view$model()              # not garbage collected


###################################################
### code chunk number 278: qt-mvc-qdf-na
###################################################
qdataFrame(model)$mpg[1] <- NA


###################################################
### code chunk number 279: qt-mvc-rtfd
###################################################
delegate <- qrTextFormattingDelegate()
view$setItemDelegate(delegate)


###################################################
### code chunk number 280: qt-mvc-table-embedded
###################################################
model <- qdataFrameModel(mtcars[,1:5])
view <- Qt$QTableView()
view$setModel(model)
window <- Qt$QWidget()
window$resize(1000, 500)
vbox <- Qt$QVBoxLayout()
vbox$addWidget(view)
window$setLayout(vbox)


###################################################
### code chunk number 281: qt-mvc-table-stretchLast
###################################################
header <- view$horizontalHeader()
header$stretchLastSection <- TRUE


###################################################
### code chunk number 282: qt-mvc-table-defaultSectionSize
###################################################
header$defaultSectionSize <- 150
header$stretchLastSection <- TRUE


###################################################
### code chunk number 283: qt-mvc-table-fit
###################################################
view$resizeColumnsToContents()


###################################################
### code chunk number 284: qt-mvc-table-stretch
###################################################
header$resizeMode(Qt$QHeaderView$Stretch)


###################################################
### code chunk number 285: qt-mvc-table-elide
###################################################
view$textElideMode <- Qt$Qt$ElideNone


###################################################
### code chunk number 286: qt-mvc-table-cascade
###################################################
header$cascadingSectionResizes <- TRUE


###################################################
### code chunk number 287: qt-mvc-lists-df
###################################################
model <- qdataFrameModel(rownames(mtcars))
view <- Qt$QListView()
view$setModel(model)


###################################################
### code chunk number 288: qt-mvc-lists-combo
###################################################
mtcars_id <- cbind(makeAndModel = rownames(mtcars), mtcars)
model <- qdataFrameModel(mtcars_id)
table_view <- Qt$QTableView()
table_view$setModel(model)
##
list_view <- Qt$QListView()
list_view$setModel(model)


###################################################
### code chunk number 289: qt-mvc-lists-sort
###################################################
DF <- qdataFrame(model)
qdataFrame(model) <- DF[order(DF$mpg),]


###################################################
### code chunk number 290: qt-mvc-lists-qslm
###################################################
model <- Qt$QStringListModel(rownames(mtcars))
list_view <- Qt$QListView()
list_view$setModel(model)


###################################################
### code chunk number 291: qt-mvc-lists-qslm-retrieve
###################################################
head(model$stringList())


###################################################
### code chunk number 292: qt-mvc-combo-box
###################################################
combo_box <- Qt$QComboBox()
combo_box$setModel(model)


###################################################
### code chunk number 293: qt-mvc-model-index-create
###################################################
index <- model$index(0, 0)
c(row = index$row(), column = index$column())


###################################################
### code chunk number 294: qt-mvc-model-index-data
###################################################
(first_car <- index$data())


###################################################
### code chunk number 295: qt-mvc-model-populate
###################################################
items <- sapply(seq(model$rowCount()), function(i) {
  model$index(i - 1, 0)$data()
})
head(items, n=6)


###################################################
### code chunk number 296: qt-mvc-model-index-data-set
###################################################
model$setData(index, toupper(first_car))


###################################################
### code chunk number 297: qt-mvc-selection-set-mode
###################################################
list_view$selectionMode <- Qt$QAbstractItemView$SingleSelection


###################################################
### code chunk number 298: qt-mvc-selection-set-behavior
###################################################
table_view$selectionBehavior <- Qt$QAbstractItemView$SelectRows


###################################################
### code chunk number 299: qt-mvc-selection-model
###################################################
selection_model <- list_view$selectionModel()


###################################################
### code chunk number 300: qt-mvc-selection-set-initial
###################################################
selection_model$select(model$index(2, 0), 
                       Qt$QItemSelectionModel$Select)


###################################################
### code chunk number 301: qt-mvc-selection-get-indices
###################################################
indices <- selection_model$selectedIndexes()
indices[[1]]$data()


###################################################
### code chunk number 302: qt-mvc-selection-multiple-set-initial-I
###################################################
sel <- Qt$QItemSelection()
sel$select(model$index(3, 0), model$index(5, 0))
sel$select(model$index(10, 0), model$index(16, 0))
sel$select(model$index(20, 0), model$index(24, 0))
selection_model$select(sel, Qt$QItemSelectionModel$ClearAndSelect)


###################################################
### code chunk number 303: qt-mvc-selection-get
###################################################
selection <- selection_model$selection()


###################################################
### code chunk number 304: qt-mvc-selection-as-list
###################################################
indicesForSelection <- function(selection) {
  selection_ranges <- as.list(selection)
  unlist(lapply(selection_ranges, function(range) {
    seq(range$top(), range$bottom())
  }))
}
indicesForSelection(selection)


###################################################
### code chunk number 305: qt-mvc-selection-changed
###################################################
selected_indices <- rep(FALSE, nrow(mtcars))
selectionChangedHandler <- function(selected, deselected) {
  selected_indices[indicesForSelection(selected)] <<- TRUE
  selected_indices[indicesForSelection(deselected)] <<- FALSE
}
qconnect(selection_model, "selectionChanged", 
         selectionChangedHandler)


###################################################
### code chunk number 306: qt-mvc-selection-set-current
###################################################
list_view$setCurrentIndex(model$index(0, 0))


###################################################
### code chunk number 307: qt-mvc-selection-select
###################################################
selection_model$select(model$index(0, 0), 
                       Qt$QItemSelectionModel$Select)


###################################################
### code chunk number 308: qt-mvc-selection-deselect
###################################################
selection_model$select(model$index(0, 0), 
                       Qt$QItemSelectionModel$Deselect)


###################################################
### code chunk number 309: qt-mvc-selection-multiple-set-initial
###################################################
selection <- Qt$QItemSelection(model$index(3, 0), 
                               model$index(10, 0))
selection_model$select(selection, 
                       Qt$QItemSelectionModel$Select)


###################################################
### code chunk number 310: qt-mvc-selection-select-row
###################################################
table_view$selectRow(1)
table_view$selectRow(2)


###################################################
### code chunk number 311: qt-mvc-selection-get-rows
###################################################
selection_model <- table_view$selectionModel()
sapply(selection_model$selectedRows(), qinvoke, "row")


###################################################
### code chunk number 312: qt-mvc-selection-select-row
###################################################
table_view$selectRow(0)


###################################################
### code chunk number 313: qt-mvc-selection-select-rows (eval = FALSE)
###################################################
## selection_model$select(selection, 
##   Qt$QItemSelectionModel$Select | Qt$QItemSelectionModel$Rows)


###################################################
### code chunk number 314: setModel
###################################################
mtcars.id <- cbind(makeAndModel = rownames(mtcars), mtcars)
model <- qdataFrameModel(mtcars.id)


###################################################
### code chunk number 315: qt-mvc-proxy
###################################################
proxy_model <- Qt$QSortFilterProxyModel()
proxy_model$setSourceModel(model)
table_view$setModel(proxy_model)
list_view$setModel(proxy_model)


###################################################
### code chunk number 316: qt-mvc-proxy-sort-enable
###################################################
table_view$sortingEnabled <- TRUE


###################################################
### code chunk number 317: qt-mvc-proxy-sort
###################################################
proxy_model$sort(1)               # mpg in column 2 of model


###################################################
### code chunk number 318: qt-mvc-proxy-filter
###################################################
proxy_model$filterKeyColumn <- 0
proxy_model$filterRegExp <- Qt$QRegExp("^Merc")


###################################################
### code chunk number 319: qt-mvc-table-hidden
###################################################
table_view$setColumnHidden(5 - 1L, TRUE)


###################################################
### code chunk number 320: qt-mvc-role-background
###################################################
mtcars_id <- cbind(makeAndModel=rownames(mtcars), mtcars)
model <- qdataFrameModel(mtcars_id, useRoles = TRUE)
qdataFrame(model)$.makeAndModel.background <- 
  list(qcolor("gray"))


###################################################
### code chunk number 321: qt-mvc-role-multi-font
###################################################
qdataFrame(model)$.mpg.hp.font <- 
  list(qfont(weight = Qt$QFont$Bold))


###################################################
### code chunk number 322: showView
###################################################
view <- Qt$QTableView()
view$setModel(model)
view$verticalHeader()$hide()       # hide default row names


###################################################
### code chunk number 323: qt-mvc-role-all-font (eval = FALSE)
###################################################
## qdataFrame(model)$.foreground <- list(qcolor("darkgray"))


###################################################
### code chunk number 324: qt-mvc-role-set-data
###################################################
list_model <- Qt$QStringListModel(rownames(mtcars))
list_model$setData(list_model$index(0, 0), "yellow", 
                   Qt$Qt$BackgroundRole)
list_view <- Qt$QListView()
list_view$setModel(list_model)


###################################################
### code chunk number 325: Widgets-MVC.Rnw:843-861
###################################################
DF <- rbind(
            c("\\code{DisplayRole}", "How data is displayed (\\class{QString})"),
            c("\\code{EditRole}", "Data for editing (\\class{QString})"),
            c("\\code{ToolTipRole}", "Displayed in tooltip (\\class{QString})"),
            c("\\code{StatusTipRole}", "Displayed in status bar (\\class{QString})"),
            c("\\code{SizeHintRole}", "Size hint for views (\\class{QSize})"),
            c("\\code{DecorationRole}", " (\\class{QColor}, \\class{QIcon}, \\class{QPixmap})"),
            c("\\code{FontRole}", "Font for default delegate (\\class{QFont})"),
            c("\\code{TextAlignmentRole}", "Alignment for default delegate (\\qtenumeration{Qt::AlignmentFlag})"),
            c("\\code{BackgroundRole}", "Background for default delegate (\\class{QBrush})"),
            c("\\code{ForegroundRole}", "Foreground for default delegate (\\class{QBrush})"),
            c("\\code{CheckStateRole}", "Indicates checked state of item (\\qtenumeration{Qt::CheckState})")
            )
colnames(DF) <- c("Constant","Description")
cat(booktabs(DF,
             colTypes=c("l","p{0.7\\textwidth}"),
             caption="Partial list of roles that an item can hold data for and the class of the data.",
             label="tab:qt:enum:itemDataRole"))


###################################################
### code chunk number 326: qt-mvc-standard-item-model
###################################################
tree_model <- Qt$QStandardItemModel(rows = 0, columns = 1)


###################################################
### code chunk number 327: qt-mvc-standard-item-set-item
###################################################
by(Cars93, Cars93$Manufacturer, function(DF) {
  tree_model$insertRow(tree_model$rowCount())
  manufacturer <- tree_model$index(tree_model$rowCount()-1L,0)
  tree_model$setData(manufacturer, DF$Manufacturer[1])
  tree_model$insertRows(0, nrow(DF), manufacturer)
  tree_model$insertColumn(0, manufacturer)
  for (i in seq_along(DF$Model)) {
    record <- tree_model$index(i-1L, 0, manufacturer)
    tree_model$setData(record, DF$Model[i])
  }
})


###################################################
### code chunk number 328: Widgets-MVC.Rnw:926-927
###################################################
tree_model <- Qt$QStandardItemModel(rows = 0, columns = 1)


###################################################
### code chunk number 329: qt-mvc-standard-item-rewrite
###################################################
by(Cars93, Cars93$Manufacturer, function(DF) {
  manufacturer <- as.character(DF$Manufacturer[1])
  manufacturer_item <- Qt$QStandardItem(manufacturer)
  tree_model$appendRow(manufacturer_item)
  children <- lapply(as.character(DF$Model), Qt$QStandardItem)
  lapply(children, manufacturer_item$appendRow)
})


###################################################
### code chunk number 330: qt-mvc-tree-view
###################################################
tree_view <- Qt$QTreeView()
tree_view$setModel(tree_model)


###################################################
### code chunk number 331: Widgets-MVC.Rnw:949-962
###################################################
## How to see model three ways
view <- Qt$QTreeView()
view$windowTitle <- "QTreeView"
view$headerHidden <- TRUE
view$setModel(tree_model); view$show(); view$raise()
##
view <- Qt$QTableView()
view$windowTitle <- "QTableView"
view$setModel(tree_model); view$show(); view$raise()
##
view <- Qt$QListView()
view$windowTitle <- "QListView"
view$setModel(tree_model); view$show(); view$raise()


###################################################
### code chunk number 332: qt-mvc-tree-view-header-hidden
###################################################
tree_view$headerHidden <- TRUE


###################################################
### code chunk number 333: ex-qt-treewidget.Rnw:9-11
###################################################
## use QTreeView to make workspace browser
require(qtbase)


###################################################
### code chunk number 334: ex-qt-treewidget.Rnw:24-114
###################################################
## From an earlier example
###################################################
### chunk number 1: ws_watcher
###################################################
library(qtbase)


###################################################
### chunk number 2: 
###################################################
qsetClass("WSWatcher", Qt$QObject, function(parent=NULL) {
  super(parent)
  updateVariables()
})


###################################################
### chunk number 3: 
###################################################
library(digest)


###################################################
### chunk number 4: 
###################################################
qsetProperty("digests", WSWatcher)


###################################################
### chunk number 5: objects_changed_property
###################################################
qsetSignal("objectsChanged", WSWatcher)


###################################################
### chunk number 6: 
###################################################
qsetProperty("objects", WSWatcher, notify="objectsChanged")


###################################################
### chunk number 7: 
###################################################
qsetProperty("old_digests", WSWatcher)
qsetProperty("old_objects", WSWatcher)


###################################################
### chunk number 8: update_variables
###################################################
qsetMethod("updateVariables", WSWatcher, function() {
  x <- sort(ls(envir=.GlobalEnv))
  objs <- sapply(mget(x, .GlobalEnv), digest)

  if((length(objs) != length(digests)) ||
     length(digests) == 0 ||
     any(objs != digests)) {
    this$old_digests <- digests         # old
    this$old_objects <- objects
    this$digests <- objs                # update cache
    this$objects <- x                   # emits signal         
  }
  invisible()
})


###################################################
### chunk number 9: change_add
###################################################
qsetMethod("changedVariables", WSWatcher, function() {
  changed <- setdiff(old_digests, digests)
  old_objects[old_digests %in% changed]
})
##
qsetMethod("addedVariables", WSWatcher, function() {
  added <- setdiff(digests, old_digests)
  objects[digests %in% added]
})





###################################################
### code chunk number 335: ex-qt-treewidget.Rnw:121-141
###################################################
addItem <- function(varname, parent_object, parent_item) {
          
  obj <- parent_object[[varname]]
  ## main interaction with tree model
  item <- Qt$QStandardItem(varname)
  class_item <- Qt$QStandardItem(paste(class(obj), 
                                      collapse = ", "))
  parent_item$appendRow(list(item, class_item))

  ## Recursively create ancestor items, if needed
  nms <- NULL
  if (is.recursive(obj)) {
    if (is.environment(obj))
      nms <- ls(obj)
    else if (!is.null(names(obj)))
      nms <- names(obj)
  }
  sapply(nms, addItem, parent_item = item, 
         parent_object = obj)
}


###################################################
### code chunk number 336: updateTopLevelItems
###################################################
updateTopLevelItems <- function(ws_watcher, view, 
                                env = .GlobalEnv) {
  ## remove these (by index)
  remove <- ws_watcher$changedVariables()
  cur_shown <- sapply(seq(model$rowCount()), 
                 function(i) model$index(i - 1, 0)$data())
  indices_to_remove <- which(cur_shown == remove)
  indices_to_remove <- sort(inds_to_remove, decreasing=TRUE)  
  ## add these (by variable name)
  new_names <- ws_watcher$addedVariables()
  
  ## replace/add these
  model <- view$model()
  view$updatesEnabled <- FALSE
  if(length(indices_to_remove))
    sapply(indices_to_remove -1L, model$removeRow)
  ## add
  sapply(new_names, addItem, parent_object = env,
         parent_item = model$invisibleRootItem())
  model$sort(0, Qt$Qt$AscendingOrder)
  view$updatesEnabled <- TRUE
}


###################################################
### code chunk number 337: ex-qt-treewidget.Rnw:181-192
###################################################
initializeTopLevelItems <- function(ws_watcher, view, 
                                    env = .GlobalEnv) 
{
   current_names <- ws_watcher$objects
   model <- view$model()
   view$updatesEnabled <- FALSE
   sapply(current_names, addItem, parent_object = env, # add
          parent_item = model$invisibleRootItem())
   model$sort(0, Qt$Qt$AscendingOrder)
   view$updatesEnabled <- TRUE
}


###################################################
### code chunk number 338: showTree
###################################################
model <- Qt$QStandardItemModel(rows = 0, columns = 2)
model$setHorizontalHeaderLabels(c("Name", "Class"))
view <- Qt$QTreeView()
view$windowTitle <- "Workspace Browser"
view$headerHidden <- FALSE
view$setModel(model)


###################################################
### code chunk number 339: initialize (eval = FALSE)
###################################################
## ws_watcher <- WSWatcher()
## ws_watcher$updateVariables()
## initializeTopLevelItems(ws_watcher, view)


###################################################
### code chunk number 340: objectsChanged (eval = FALSE)
###################################################
## qconnect(ws_watcher, "objectsChanged", function() 
##          updateTopLevelItems(ws_watcher, view))


###################################################
### code chunk number 341: taskCallback (eval = FALSE)
###################################################
## ## add callback
## id <- addTaskCallback(function(expr, value, ok, visible) {
##   ws_watcher$updateVariables()
##   TRUE
## })
## ## view
## view$show()
## view$raise()


###################################################
### code chunk number 342: qt-mvc-check-editable
###################################################
(tree_model$index(0, 0)$flags() & Qt$Qt$ItemIsEditable) > 0


###################################################
### code chunk number 343: qt-mvc-edit-analyze
###################################################
DF <- mtcars 
DF$Analyze <- TRUE
model <- qdataFrameModel(DF, editable = "Analyze")


###################################################
### code chunk number 344: qt-mvc-edit-trigger
###################################################
view$editTriggers <- Qt$QAbstractItemView$NoEditTriggers


###################################################
### code chunk number 345: qt-mvc-dnd-dragenabled
###################################################
view$dragEnabled <- TRUE


###################################################
### code chunk number 346: qt-mvc-dnd-acceptdrops
###################################################
view$acceptDrops <- TRUE
view$showDropIndicator <- TRUE


###################################################
### code chunk number 347: qt-mvc-dnd-internalmove
###################################################
view$dragDropMode <- Qt$QAbstractItemView$InternalMove


###################################################
### code chunk number 348: ex-qt-dnd-table.Rnw:1-3
###################################################
library(qtbase)
#rm(list=ls())                           # clear out


###################################################
### code chunk number 349: qt-mvc-dnd-variable-selector
###################################################
qsetClass("VariableSelector", Qt$QWidget, 
          function(parent = NULL) {
  super(parent)
  ## widgets
  this$df_combo_box <- Qt$QComboBox()
  this$variable_list <- Qt$QListView()
  this$variable_list$setModel(
       qdataFrameModel(data.frame(), this, 
                       useRoles = TRUE))
  this$variable_list$dragEnabled <- TRUE

  ## layout
  layout <- Qt$QVBoxLayout()
  layout$addWidget(df_combo_box)
  layout$addWidget(variable_list)
  variable_list$setSizePolicy(Qt$QSizePolicy$Expanding, 
                        Qt$QSizePolicy$Expanding)
  setLayout(layout)
  
  updateDataSets()
  qconnect(df_combo_box, "activated(int)", function(ind) {
    this$dataFrame <- df_combo_box$currentText
  })
})


###################################################
### code chunk number 350: qt-mvc-dnd-update-datasets
###################################################
qsetMethod("updateDataSets", VariableSelector, function() {
  current_text <- df_combo_box$currentText
  df_combo_box$clear()
  DFs <- ProgGUIinR:::avail_dfs(.GlobalEnv)
  if(length(DFs)) {
    this$df_combo_box$addItems(DFs)
    if(is.null(current_text) || !current_text %in% DFs) {
      this$df_combo_box$currentIndex <- -1
      this$dataFrame <- NULL
    } else {
      this$df_combo_box$currentIndex <- 
        which(current_text == DFs)
      this$dataFrame <- current_text
    }
  }
})


###################################################
### code chunk number 351: ex-qt-dnd-table.Rnw:72-83
###################################################
## The \function{getIcon} helper function provides an icon from the class of a
## column:
require(grid)
getIcon <- function(x) {
  f <- tempfile()
  png(file=f, width=16, height=16)
  grid::grid.newpage()
  grid::grid.draw(ProgGUIinR:::make_icon(x))
  dev.off()
  Qt$QIcon(f)
}


###################################################
### code chunk number 352: qt-mvc-dnd-dataset
###################################################
qsetProperty("dataFrame", VariableSelector, 
             write = function(DF) {
               if (is.null(DF))
                 DF <- data.frame()
               else if (is.character(DF)) 
                 DF <- get(DF, .GlobalEnv)
               ##
               model <- variable_list$model()
               icons <- lapply(DF, getIcon)
               qdataFrame(model) <- 
                 data.frame(variable=names(DF),
                            variable.decoration=I(icons))
               this$.dataFrame <- DF
               dataFrameChanged()
             })


###################################################
### code chunk number 353: qt-mvc-dnd-datasetChanged
###################################################
qsetSignal("dataFrameChanged", VariableSelector)


###################################################
### code chunk number 354: DropLabelRotation
###################################################
qsetClass("VariableLabel", Qt$QLabel, function(parent=NULL) {
  super(parent)
  this$rotation <- 0L
  setAcceptDrops(TRUE)
  setAlignment(Qt$Qt$AlignHCenter | Qt$Qt$AlignVCenter)
})


###################################################
### code chunk number 355: qt-mvc-dnd-rotation
###################################################
qsetProperty("rotation", VariableLabel)
qsetProperty("variable_name", VariableLabel)


###################################################
### code chunk number 356: qt-mvc-dnd-drop
###################################################
qsetSignal("variableNameDropped", VariableLabel)


###################################################
### code chunk number 357: qt-mvc-dnd-get-variable-name
###################################################
variableNameFromMimeData <- function(mime_data) {
  name <- NULL
  RDA_MIME_TYPE <- "application/x-rlang-transport"
  if(mime_data$hasFormat(RDA_MIME_TYPE)) {
    name_list <- unserialize(mime_data$data(RDA_MIME_TYPE))
    if (length(name_list) && is.character(name_list[[1]]))
      name <- name_list[[1]]
  }
  name
}


###################################################
### code chunk number 358: ex-qt-dnd-table.Rnw:159-170
###################################################
qsetMethod("dragEnterEvent", VariableLabel, function(event) {
  mime_data <- event$mimeData()
  if(!is.null(variableNameFromMimeData(mime_data))) {
    setForegroundRole(Qt$QPalette$Dark)
    event$acceptProposedAction()
  }
})
qsetMethod("dragLeaveEvent", VariableLabel, function(event) {
  setForegroundRole(Qt$QPalette$WindowText)
  event$accept()
})


###################################################
### code chunk number 359: dropEvent
###################################################
qsetMethod("dropEvent", VariableLabel, function(event) {
  setForegroundRole(Qt$QPalette$WindowText)  
  mime_data <- event$mimeData()
  this$variable_name <- variableNameFromMimeData(mime_data)
  if(!is.null(variable_name)) {
    this$text <- variable_name
    variableNameDropped()
    setBackgroundRole(Qt$QPalette$Window)
    event$acceptProposedAction()
  }
})


###################################################
### code chunk number 360: ex-qt-dnd-table.Rnw:195-208
###################################################
qsetMethod("paintEvent", VariableLabel, function(event) {
  painter <- Qt$QPainter()
  painter$begin(this)
  
  painter$save()
  painter$translate(width / 2, height / 2)
  painter$rotate(-(rotation))
  rect <- painter$boundingRect(0, 0, 0, 0, 
                               Qt$Qt$AlignCenter, text)
  painter$drawText(rect, Qt$Qt$AlignCenter, text)
  painter$restore()
  painter$end()
})


###################################################
### code chunk number 361: XtabsWidget
###################################################
qsetClass("XtabsWidget", Qt$QWidget, function(parent = NULL) {
  super(parent)
  initWidgets()
  initLayout()
})


###################################################
### code chunk number 362: initWidgets
###################################################
qsetMethod("initWidgets", XtabsWidget, function() {
  this$xlabel <- VariableLabel()
  qconnect(xlabel, "variableNameDropped", invokeXtabs)

  this$ylabel <- VariableLabel()
  pt <- ylabel$font$pointSize()
  ylabel$minimumWidth <- 2*pt; ylabel$maximumWidth <- 2*pt
  ylabel$rotation <- 90L
  qconnect(ylabel, "variableNameDropped", invokeXtabs)
  
  this$table_view <- Qt$QTableView()
  table_view$setModel(qdataFrameModel(data.frame(), this))
  clearLabels()
})


###################################################
### code chunk number 363: ex-qt-dnd-table.Rnw:252-261
###################################################
## Not shown
qsetMethod("clearLabels", XtabsWidget, function() {
  point_size <- xlabel$font$pointSize()
  xlabel$text <- "Drop x factor here"
  xlabel$minimumHeight <- 2 * point_size
  
  ylabel$text <- "Drop y factor here"
  ylabel$minimumWidth <- 2 * point_size
})


###################################################
### code chunk number 364: ex-qt-dnd-table.Rnw:264-275
###################################################
## Not shown
qsetMethod("initLayout", XtabsWidget, function() {
  layout <- Qt$QGridLayout()
  setLayout(layout)
  layout$addWidget(xlabel, 0, 1, 1, 3)
  layout$addWidget(ylabel, 1, 0, 1, 1)
  layout$addWidget(table_view, 1, 1, 1, 3)
  
  layout$setColumnStretch(2, 1)
  layout$setRowStretch(1, 1)
})


###################################################
### code chunk number 365: ex-qt-dnd-table.Rnw:280-296
###################################################
## Hide call to xtabs
## Return NULL if not okay, otherwise a table object
call_xtabs <- function(DF, x, y) {
  if(is.character(DF))
    DF <- get(DF)
  if(is.null(x)) {
    table <- NULL
  } else if(is.null(y)) {
    f <- formula(sprintf("~ %s", x))
    table <- xtabs(f, data = DF)
  } else { 
    f <- formula(sprintf("~ %s + %s", y, x))
    table <- xtabs(f, data = DF)
  } 
  table
}


###################################################
### code chunk number 366: makeTable
###################################################
qsetMethod("invokeXtabs", XtabsWidget, function() {
  if (is.null(dataFrame))
    return()

  x <- xlabel$variable_name
  y <- ylabel$variable_name
  
  if(!is.null(table <- call_xtabs(dataFrame, x, y)))
     updateTableView(table)
})


###################################################
### code chunk number 367: updateTableWidget
###################################################
qsetMethod("updateTableView", XtabsWidget, function(table) {
  model <- table_view$model()
  if (length(dim(table)) == 1)
    qdataFrame(model) <- data.frame(count = unclass(table))
  else qdataFrame(model) <- data.frame(unclass(table))
})


###################################################
### code chunk number 368: qt-mvc-dnd-dataframe-xtabs
###################################################
qsetProperty("dataFrame", XtabsWidget, 
             write = function(dataFrame) { 
               clearLabels()
               this$.dataFrame <- dataFrame
             })


###################################################
### code chunk number 369: ex-qt-dnd-table.Rnw:335-337
###################################################
## Not shown
require(MASS); data(Cars93); data(Aids2)


###################################################
### code chunk number 370: ex-qt-dnd-table.Rnw:339-348
###################################################
w <- Qt$QSplitter()
w$setWindowTitle("GUI for xtabs()")
w$addWidget(vs <- VariableSelector())
w$addWidget(tw <- XtabsWidget())
w$setStretchFactor(1, 1)
qconnect(vs, "dataFrameChanged", function() {
  tw$dataFrame <- vs$dataFrame
})
w$show(); w$raise()


###################################################
### code chunk number 371: qt-mvc-listwidget-additems
###################################################
list_widget <- Qt$QListWidget()
list_widget$addItems(state.name)


###################################################
### code chunk number 372: qt-mvc-listwidget-item
###################################################
item <- Qt$QListWidgetItem("Puerto Rico", list_widget)


###################################################
### code chunk number 373: qt-mvc-listwidget-itemat
###################################################
first <- list_widget$item(0)
first$text()


###################################################
### code chunk number 374: qt-mvc-listwidget-selectionmode
###################################################
list_widget$selectionMode <- Qt$QListWidget$ExtendedSelection


###################################################
### code chunk number 375: qt-mvc-listwidget-select
###################################################
sapply(grep("^A", state.name), 
       function(i) list_widget$item(i - 1)$setSelected(TRUE))


###################################################
### code chunk number 376: qt-mvc-listwidget-selected
###################################################
selected_items <- list_widget$selectedItems()
sapply(selected_items, qinvoke, "text")


###################################################
### code chunk number 377: qt-mvc-listwidget-selectionchanged
###################################################
qconnect(list_widget, "itemSelectionChanged", function() {
  selected <- list_widget$selectedItems()
  selected_text <- sapply(selected, qinvoke, "text")
  message("Selected: ", paste(selected_text, collapse = ", "))
})


###################################################
### code chunk number 378: qt-mvc-listwidget-checked
###################################################
items <- sapply(seq(list_widget$count) - 1L, list_widget$item)
sapply(items, qinvoke, "setCheckState", Qt$Qt$Unchecked)
## check selected
selected <- list_widget$selectedItems()
sapply(selected, function(x) x$setCheckState(Qt$Qt$Checked))
## clear selection now
list_widget$selectionModel()$clear()
list_widget$selectionMode <- Qt$QListWidget$NoSelection


###################################################
### code chunk number 379: Widgets-MVC.Rnw:1224-1226
###################################################
state <- sapply(items, "qinvoke", "checkState")
head(state, n = 8)                     # 2 is checked, 0 not


###################################################
### code chunk number 380: ex-qt-edit-data-frame.Rnw:20-28
###################################################
## Create a model for displaying a data frame -- like qdataFrameModel -- only
## with more Qt-like control over the display of columnms
## The key methods are data(), setData() and flags()
##
## For performance reasons -- which are considerable here and the main reason
## we discourage this approach -- 
## we avoid header Data and put variable names in first row!
library(qtbase)


###################################################
### code chunk number 381: DfModel
###################################################
qsetClass("DfModel", Qt$QAbstractTableModel, 
          function(DF = data.frame(V1 = character(0)), 
                   parent = NULL) 
          {
            super(parent)
            this$DF <- DF
          })


###################################################
### code chunk number 382: ex-qt-edit-data-frame.Rnw:46-50
###################################################
qsetProperty("DF", DfModel, write = function(DF) {
  this$.DF <- DF
  dataChanged(index(0, 0), index(nrow(DF), ncol(DF)))
})


###################################################
### code chunk number 383: ex-qt-edit-data-frame.Rnw:56-60
###################################################
qsetMethod("rowCount", DfModel, 
           function(index) nrow(this$DF) + 1)
qsetMethod("columnCount", DfModel, 
           function(index) ncol(this$DF))


###################################################
### code chunk number 384: displayRoleMethod
###################################################
display_role <- function(x, row, ...) UseMethod("display_role")
display_role.default <- function(x, row) 
  sprintf("%s", x[row])
display_role.numeric <- function(x, row) 
  sprintf("%.2f", x[row])
display_role.integer <- function(x, row) 
  sprintf("%d", x[row])


###################################################
### code chunk number 385: ex-qt-edit-data-frame.Rnw:86-125
###################################################
## Not shown
## text alignment to indicate to user different types of numeric values
textAlignmentRole <- function(x, row, context) UseMethod("textAlignmentRole")
textAlignmentRole.default <- function(x, row, context) Qt$Qt$AlignCenter
textAlignmentRole.integer <- function(x, row, context) Qt$Qt$AlignRight
textAlignmentRole.numeric <- function(x, row, context) Qt$Qt$AlignRight

## sets the background color
## returns a QBrush object
backgroundRole <- function(x, row, context) UseMethod("backgroundRole")
backgroundRole.default <- function(x, row, context) Qt$QBrush(Qt$QColor(0,0,0,0)) # transparent
backgroundRole.factor <- function(x, row, context) Qt$QBrush(Qt$QColor("yellow"))

##' Size hint role
##' XXX Doesn't get called
sizeHintRole <- function(x, row, context) UseMethod("sizeHintRole")
sizeHintRole.default <- function(x, row, context) {
  sz <- max(sapply(x, function(i) nchar(as.character(i))))
  avg <- Qt$QFontMetric(Qt$QFont())$averageCharWidth()
  Qt$QSize(sz*avg, sample(c(25L, 50L), 1))
}

## some text for a tooltip
toolTipRole <- function(x, row, context) UseMethod("toolTipRole")
toolTipRole.default <- function(x, row, context)
  sprintf("Value in vector with class %s", paste(class(x), collapse=","))
toolTipRole.factor <- function(x, row, context) {
  x <- levels(x)
  n <- length(x)
  p <- 6
  q <- n %/% p
  r <- n %% p
  l <- paste(paste(apply(matrix(x[1:(p*q)], byrow=T, ncol=6), 1, paste, collapse=","), collapse=",\n"),
        paste(x[(n-r+1):n], collapse=","), sep=",\n")
  sprintf("Factor with levels:\n%s", l)
}
toolTipRole.logical <- function(x, row, context) sprintf("Logical vector")

whatsThisRole <- toolTipRole


###################################################
### code chunk number 386: ex-qt-edit-data-frame.Rnw:130-148 (eval = FALSE)
###################################################
## qsetMethod("data", DfModel, function(index, role) {
##   row <- index$row()
##   col <- index$column() + 1
## 
##   if(role == Qt$Qt$DisplayRole) {
##     if(row > 0)
##       display_role(DF[,col], row)
##     else
##       names(DF)[col]
##   } else if(role == Qt$Qt$EditRole) {
##     if(row > 0)
##       as.character(DF[row, col])
##     else
##       names(DF)[col]
##   } else {
##     NULL
##   }
## })


###################################################
### code chunk number 387: data
###################################################
## this is not shown in text, but is the definition of the data method
qsetMethod("data", DfModel, function(index, role) {
  row <- index$row()
  col <- index$column() + 1

  
  if(role == Qt$Qt$DisplayRole) {
    if(row > 0)
      display_role(DF[,col], row)
    else
      names(DF)[col]
  } else if(role == Qt$Qt$EditRole) {
    if(row > 0)
      as.character(DF[row, col])
    else
      names(DF)[col]
  } else if(role == Qt$Qt$TextAlignmentRole) {
    if(row > 0)
      textAlignmentRole(DF[, col], row)
    else
      Qt$Qt$AlignCenter | Qt$Qt$AlignTop
  } else if(role == Qt$Qt$BackgroundRole) {
    if(row > 0)
      backgroundRole(DF[, col], row)
    else
      Qt$QBrush(Qt$QColor("gray"))
  } else if(role == Qt$Qt$ForegroundRole) {
    if(row == 0)
      Qt$QBrush(Qt$QColor("white"))
  } else if(role == Qt$Qt$ToolTipRole) {
    if(row > 0)
      toolTipRole(DF[,col], row)
    else
      ""
  } else if(role == Qt$Qt$WhatsThisRole) {
    if(row > 0)
      whatsThisRole(DF[,col], row)
    else
      ""
  } else if(role == Qt$Qt$SizeHintRole) {
    if(row > 0)
      sizeHintRole(DF[,col], row)
    else
      NULL
  } else {
    NULL
  }
})


###################################################
### code chunk number 388: ex-qt-edit-data-frame.Rnw:205-213
###################################################
qsetMethod("flags", DfModel, function(index) {
  if(!index$isValid()) {
    return(Qt$Qt$ItemIsEnabled)
  } else {
    current_flags <- super("flags", index)
    return(current_flags | Qt$Qt$ItemIsEditable)
  }
})


###################################################
### code chunk number 389: fit_in
###################################################
fit_in <- function(x, value) UseMethod("fit_in")
fit_in.default <- function(x, value) value
fit_in.numeric <- function(x, value) as.numeric(value)


###################################################
### code chunk number 390: notShown
###################################################
## more methods for fit in
fit_in.integer <- function(x, value) as.integer(value)
fit_in.logical <- function(x, value) {
  if(toupper(value) %in% c("T","F","TRUE","FALSE")) {
    as.logical(value)
  } else {
    as.logical(as.numeric(value))
  }
}


###################################################
### code chunk number 391: ex-qt-edit-data-frame.Rnw:241-261
###################################################
qsetMethod("setData", DfModel, function(index, value, role) {
  if(index$isValid() && role == Qt$Qt$EditRole) {
    DF <- this$DF
    row <- index$row()
    col <- index$column() + 1

    if(row > 0) {
      x <- DF[, col]
      DF[row, col] <- fit_in(x, value)
    } else {
      names(DF)[col] <- value
    }
    this$DF <- DF
    dataChanged(index, index)

    return(TRUE)
  } else {
     super("setData", index, value, role)
  }
})


###################################################
### code chunk number 392: ex-qt-edit-data-frame.Rnw:267-280
###################################################
qsetMethod("setColumn", DfModel, function(col, value) {
  ## pad with NA if needed
  n <- nrow(this$DF)
  if(length(value) < n)
    value <- c(value, rep(NA, n - length(value)))
  value <- value[1:n]
  DF <- this$DF
  DF[,col] <- value
  this$DF <- DF         # only notify about this column
  dataChanged(index(0, col - 1), 
              index(rowCount() - 1, col - 1))
  return(TRUE)
})


###################################################
### code chunk number 393: addColumn
###################################################
qsetMethod("addColumn", DfModel, function(name, value) {
  DF <- this$DF
  if(name %in% names(DF)) {
    return(setColumn(min(which(name == names(DF))), value))
  }  
  beginInsertColumns(Qt$QModelIndex(),
                     columnCount(), columnCount())
  DF[[name]] <- value
  this$DF <- DF
  endInsertColumns()
  return(TRUE)
})


###################################################
### code chunk number 394: ex-qt-edit-data-frame.Rnw:302-306
###################################################
model <- DfModel(mtcars)

view <- Qt$QTableView()
view$setModel(model)


###################################################
### code chunk number 395: customizeView
###################################################
trigger_flag <- Qt$QAbstractItemView$DoubleClicked | 
                Qt$QAbstractItemView$SelectedClicked |
                Qt$QAbstractItemView$EditKeyPressed
view$setEditTriggers(trigger_flag)
view$verticalHeader()$setHidden(TRUE)
view$horizontalHeader()$setHidden(TRUE)


###################################################
### code chunk number 396: ex-qt-edit-data-frame.Rnw:319-321
###################################################
view$show()
view$raise()


###################################################
### code chunk number 397: qt-mvc-mapper-map
###################################################
data(Cars93, package="MASS")
model <- qdataFrameModel(Cars93, editable=names(Cars93))
mapper <- Qt$QDataWidgetMapper()
mapper$setModel(model)
##
label <- Qt$QLineEdit()
mapper$addMapping(label, 1)


###################################################
### code chunk number 398: qt-mvc-mapper-select
###################################################
table_view <- Qt$QTableView()
table_view$setModel(model)
qconnect(table_view$selectionModel(), "currentRowChanged", 
         function(cur,prev) mapper$setCurrentIndex(cur$row()))


###################################################
### code chunk number 399: qt-mvc-mapper-layout
###################################################
window <- Qt$QWidget()
layout <- Qt$QVBoxLayout()
window$setLayout(layout)
layout$addWidget(table_view)
layout$addWidget(label)


###################################################
### code chunk number 400: ex-qt-custom-view.Rnw:17-18
###################################################
library(qtbase)


###################################################
### code chunk number 401: CustomView
###################################################
qsetClass("MeanLabel", Qt$QLabel, 
          function(model, column = 0, parent = NULL) 
          {
            super(parent)
            this$model <- model
            this$column <- column
            updateMean()     # initialize text
            qconnect(model, "dataChanged", 
                     function(top_left, bottom_right) {
                       if (top_left$column() <= column && 
                           bottom_right$column() >= column)
                         updateMean()
                     })
          })


###################################################
### code chunk number 402: label
###################################################
qsetMethod("updateMean", MeanLabel, function() {
  if(is.null(model)) {
    text <- "No model"
  } else {
    DF <- qdataFrame(model)
    colname <- colnames(DF)[column + 1L]
    text <- sprintf("Mean for '%s': %s", colname, 
                    mean(DF[,colname]))
  }
  this$text <- text
}, access="private")


###################################################
### code chunk number 403: testItOut
###################################################
model <- qdataFrameModel(mtcars, editable = colnames(mtcars))

table_view <- Qt$QTableView()
table_view$setModel(model)
table_view$setEditTriggers(Qt$QAbstractItemView$DoubleClicked)
##
mean_label <- MeanLabel(model)
##
window <- Qt$QWidget()
layout <- Qt$QVBoxLayout()
window$setLayout(layout)
layout$addWidget(table_view)
layout$addWidget(mean_label)


###################################################
### code chunk number 404: ex-qt-custom-view.Rnw:74-76
###################################################
window$show()
window$raise()


###################################################
### code chunk number 405: QTextEdit
###################################################
text_edit <- Qt$QTextEdit()


###################################################
### code chunk number 406: Widgets-MVC.Rnw:1350-1352 (eval = FALSE)
###################################################
## ## not shown
## text_edit$show(); text_edit$raise()


###################################################
### code chunk number 407: Widgets-MVC.Rnw:1363-1365
###################################################
text_edit$setPlainText("The quick brown fox")
text_edit$append("jumped over the lazy dog")


###################################################
### code chunk number 408: Widgets-MVC.Rnw:1376-1377
###################################################
text_edit$toPlainText()


###################################################
### code chunk number 409: qt-mvc-textedit-cursor
###################################################
n <- nchar(text_edit$toPlainText())
cursor <- text_edit$textCursor()
cursor$setPosition(n)
text_edit$setTextCursor(cursor)


###################################################
### code chunk number 410: qt-mvc-textedit-image
###################################################
cursor$setPosition(0)                   # move to beginning
style <- Qt$QApplication$style()
icon <- style$standardIcon(Qt$QStyle$SP_DialogOkButton)
sz <- qsize(32L,32L)
anImage <- icon$pixmap(icon$actualSize(sz))$toImage()
cursor$insertImage(anImage)


###################################################
### code chunk number 411: Widgets-MVC.Rnw:1443-1447
###################################################
cursor <- text_edit$textCursor()
cursor$movePosition(Qt$QTextCursor$Start) # MoveAnchor default
cursor$movePosition(Qt$QTextCursor$Down)  # down one line
text_edit$setTextCursor(cursor)


###################################################
### code chunk number 412: qt-mvc-textedit-selection
###################################################
text_edit$textCursor()$selectedText()  # no current selection


###################################################
### code chunk number 413: Widgets-MVC.Rnw:1474-1480
###################################################
cursor <- Qt$QTextCursor(text_edit$document()) 
cursor$movePosition(Qt$QTextCursor$Start)     # as before
cursor$movePosition(Qt$QTextCursor$Down)      # moves anchor
cursor$movePosition(Qt$QTextCursor$WordRight, # anchor fixed
                    Qt$QTextCursor$KeepAnchor, 3)
text_edit$setTextCursor(cursor)


###################################################
### code chunk number 414: qt-mvc-textedit-selection
###################################################
cursor$selectedText()


###################################################
### code chunk number 415: qt-mvc-textedit-selectionChanged
###################################################
qconnect(text_edit, "textChanged", function() {
  message("Text has changed to", text_edit$toPlainText())
})
##
qconnect(text_edit, "cursorPositionChanged", function() {
  message("Cursor has changed. It is now in position", 
          text_edit$textCursor()$position())
})
##
qconnect(text_edit, "selectionChanged", function() {
  message("text: ", text_edit$textCursor()$selectedText())
})


###################################################
### code chunk number 416: qt-mvc-textedit-wrap
###################################################
text_edit$lineWrapMode <- Qt$QTextEdit$NoWrap


###################################################
### code chunk number 417: LoremIpsum
###################################################
## in misc.R
LoremIpsum <- paste(
"Lorem ipsum dolor sit amet, consectetur adipisicing elit,",
"sed do eiusmod tempor incididunt ut labore et dolore magna aliqua.",
"Ut enim ad minim veniam, quis nostrud exercitation ullamco laboris nisi",
"ut aliquip ex ea commodo consequat. Duis aute irure dolor in reprehenderit",
"in voluptate velit esse cillum dolore eu fugiat nulla pariatur. Excepteur",
"sint occaecat cupidatat non proident, sunt in culpa qui officia deserunt mollit",
"anim id est laborum.",
"\n",
"Sed ut perspiciatis unde omnis iste natus error sit voluptatem accusantium doloremque",
"laudantium, totam rem aperiam, eaque ipsa quae ab illo inventore veritatis et quasi",
"architecto beatae vitae dicta sunt explicabo. Nemo enim ipsam voluptatem quia voluptas",
"sit aspernatur aut odit aut fugit, sed quia consequuntur magni dolores eos qui ratione",
"voluptatem sequi nesciunt. Neque porro quisquam est, qui dolorem ipsum quia dolor sit amet,",
"consectetur, adipisci velit, sed quia non numquam eius modi tempora incidunt ut labore et",
"dolore magnam aliquam quaerat voluptatem. Ut enim ad minima veniam, quis nostrum",
"exercitationem ullam corporis suscipit laboriosam, nisi ut aliquid ex ea commodi",
"consequatur? Quis autem vel eum iure reprehenderit qui in ea voluptate velit esse quam",
"nihil molestiae consequatur, vel illum qui dolorem eum fugiat quo voluptas nulla pariatur?",
sep="\n")


###################################################
### code chunk number 418: findExample
###################################################
text_edit <- Qt$QTextEdit(LoremIpsum)          # some text
text_edit$find("qui", Qt$QTextDocument$FindWholeWords)
text_edit$textCursor()$selection()$toPlainText()


###################################################
### code chunk number 419: QTextEditWithMenu
###################################################
qsetClass("QTextEditWithCompletions", Qt$QTextEdit)
##
qsetMethod("contextMenuEvent", QTextEditWithCompletions, 
           function(event) 
           {
  menu <- this$createStandardContextMenu()
  if(this$textCursor()$hasSelection()) {
    selection <- this$textCursor()$selectedText()
    completions <- utils:::matchAvailableTopics(selection)
    completions <- setdiff(completions, selection)
    if(length(completions) > 0 && length(completions) < 25) {
      menu$addSeparator()                  # add actions
      sapply(completions, function(completion) {
        action <- Qt$QAction(completion, this)
        qconnect(action, "triggered", function(checked) {
          insertPlainText(completion)
        })
        menu$addAction(action)
      })
    }
  }
  menu$exec(event$globalPos())
})
text_edit <- QTextEditWithCompletions()


###################################################
### code chunk number 420: raise (eval = FALSE)
###################################################
## text_edit$show()
## text_edit$raise()


###################################################
### code chunk number 421: Widgets-MVC.Rnw:1655-1706 (eval = FALSE)
###################################################
## ## Not shown
## 
## ##' start with command, return output in HTML
## runCommand <- function(chunk) {
##   require(R2HTML)
##   f <- tempfile()
##   chunkexps <- try(parse(text=chunk), silent=TRUE)
##   if(inherits(chunkexps, "try-error")) 
##   stop("Error")
##   
##   out <- ""
##   for(i in seq_along(chunkexps)) {
##     chunk <- chunkexps[[i]]
##     out <- paste(out, "<h3>", paste(capture.output(chunk), collapse="<br>"), "</h3><br />",  sep="")
##     HTML(eval(chunk), file=f, append=FALSE)
##     out <- paste(out, paste(readLines(f, warn=FALSE), collapse="<br>"), sep="")
##   }
##   out
## }
## 
## 
## 
## 
## w <- Qt$QGroupBox("Simple CLI")
## lyt <- Qt$QGridLayout()
## w$setLayout(lyt)
## 
## cli <- Qt$QLineEdit()
## out <- Qt$QTextEdit()
## 
## 
## lyt$addWidget(Qt$QLabel("Command:"), 0, 0)
## lyt$addWidget(cli, 0, 1)
## lyt$addWidget(Qt$QLabel("Output:"), 1, 0, Qt$Qt$AlignTop)
## lyt$addWidget(out, 1, 1)
## 
## lyt$setRowStretch(1,1)
## lyt$setColumnStretch(1,1)
## 
## qconnect(cli, "editingFinished", function() {
##   chunk <- cli$text
##   htmlized <- runCommand(chunk)
##   out$setHtml(htmlized)
##   cli$setFocus(Qt$Qt$ActiveWindowFocusReason)
##   cli$setSelection(0, nchar(cli$text))
## })
## 
## w$setMinimumSize(800,400)
## w$show()
## w$raise()
## 


###################################################
### code chunk number 422: qt-app-construct
###################################################
main_window <- Qt$QMainWindow()


###################################################
### code chunk number 423: qt-app-table
###################################################
data(mtcars)
model <- qdataFrameModel(mtcars, editable = TRUE)
table_view <- Qt$QTableView()
table_view$setModel(model)
main_window$setCentralWidget(table_view)


###################################################
### code chunk number 424: qt-app-action-construct
###################################################
open_action <- Qt$QAction("Open", main_window)


###################################################
### code chunk number 425: qt-app-action-status-tip
###################################################
open_action$statusTip <- "Load a spreadsheet from a CSV file"


###################################################
### code chunk number 426: qt-app-action-icon
###################################################
style <- Qt$QApplication$style()
open_action$icon <- 
  style$standardIcon(Qt$QStyle$SP_DialogOpenButton)


###################################################
### code chunk number 427: qt-app-action-trigger
###################################################
qconnect(open_action, "triggered", function() {
  filename <- Qt$QFileDialog$getOpenFileName()
  table_view$model <- 
    qdataFrameModel(read.csv(filename), editable = TRUE)
})


###################################################
### code chunk number 428: qt-app-action-toggle
###################################################
save_on_exit_action <- Qt$QAction("Save on exit", main_window)
save_on_exit_action$checkable <- TRUE


###################################################
### code chunk number 429: qt-app-action-radio
###################################################
just_group <- Qt$QActionGroup(main_window)
just_action <- list()
just_action$left <- Qt$QAction("Left Align", just_group)
just_action$right <- Qt$QAction("Right Align", just_group)
just_action$center <- Qt$QAction("Center", just_group)
sapply(just_action, function(action) action$checkable <- TRUE)


###################################################
### code chunk number 430: ApplicationWindow.Rnw:130-137
###################################################
sapply(just_action, function(action) 
       qconnect(action, "changed", function() {
         button_number <- 
           which(sapply(just_action, `[[`, "checked"))
         message("Button ", button_number, " was depressed")
       })
       )


###################################################
### code chunk number 431: actionGroupSignal
###################################################
qconnect(just_group, "triggered", function(action) {
  message(action$text)
})


###################################################
### code chunk number 432: qt-app-action-shortcut
###################################################
open_action$setShortcut(Qt$QKeySequence(Qt$QKeySequence$Open))


###################################################
### code chunk number 433: qt-app-menubar
###################################################
menubar <- Qt$QMenuBar()
main_window$setMenuBar(menubar)


###################################################
### code chunk number 434: qt-app-menubar-addMenu
###################################################
file_menu <- Qt$QMenu("File")
menubar$addMenu(file_menu)
edit_menu <- Qt$QMenu("Edit")
menubar$addMenu(edit_menu)


###################################################
### code chunk number 435: qt-app-menubar-populate
###################################################
file_menu$addAction(open_action)
file_menu$addSeparator()
file_menu$addAction(save_on_exit_action)
file_menu$addSeparator()
quit_action <- file_menu$addAction("Quit")

just_menu <- edit_menu$addMenu("Justification")
sapply(just_action, just_menu$addAction)


###################################################
### code chunk number 436: qt-app-context-add-widgets
###################################################
sort_menu <- Qt$QMenu("Sort by")
sapply(colnames(qdataFrame(model)), sort_menu$addAction)
table_view$addAction(sort_menu$menuAction())


###################################################
### code chunk number 437: qt-app-context-policy
###################################################
table_view$contextMenuPolicy <- Qt$Qt$ActionsContextMenu


###################################################
### code chunk number 438: qt-app-context-signal
###################################################
showCompletionPopup <- function(event, edit) {
  popup <- Qt$QMenu()
  completions <- utils:::matchAvailableTopics(ed$text)
  completions <- head(completions, 10) # trim if large
  sapply(completions, function(completion) {
    action <- popup$addAction(completion)
    qconnect(action, "triggered", 
             function(...) edit$setText(completion))
  })
  popup$popup(edit$mapToGlobal(qpoint(0L,0L)))
}
##
edit <- Qt$QLineEdit()
edit$contextMenuPolicy <- Qt$Qt$CustomContextMenu
qconnect(edit, "customContextMenuRequested", 
         showCompletionPopup, user.data = edit)


###################################################
### code chunk number 439: qt-app-toolbar
###################################################
toolbar <- Qt$QToolBar()
main_window$addToolBar(toolbar)


###################################################
### code chunk number 440: ApplicationWindow.Rnw:336-343
###################################################
## Before adding some actions to our toolbar, we define a function
## \code{getIcon} that loads a \class{QIcon} from a file in the
## \pkg{gWidgets} package:
getIcon <- function(nm) {
  fname <- system.file(sprintf("images/%s.gif", nm), package="gWidgets")
  Qt$QIcon(fname)
}


###################################################
### code chunk number 441: ApplicationWindow.Rnw:348-359
###################################################
file_actions <- list()
file_actions$open <- Qt$QAction("Open", main_window)
file_actions$open$setIcon(getIcon("open"))
file_actions$save <- Qt$QAction("Save", main_window)
file_actions$save$setIcon(getIcon("save"))

plot_actions <- list()
plot_actions$barplot <- Qt$QAction("Barplot", main_window)
plot_actions$barplot$setIcon(getIcon("barplot"))
plot_actions$boxplot <- Qt$QAction("Boxplot", main_window)
plot_actions$boxplot$setIcon(getIcon("boxplot"))


###################################################
### code chunk number 442: ApplicationWindow.Rnw:364-367
###################################################
sapply(file_actions, toolbar$addAction)
toolbar$addSeparator()
sapply(plot_actions, toolbar$addAction)


###################################################
### code chunk number 443: ApplicationWindow.Rnw:369-371
###################################################
main_window$show()
main_window$raise()


###################################################
### code chunk number 444: qt-app-toolbar-style
###################################################
toolbar$setToolButtonStyle(Qt$Qt$ToolButtonTextUnderIcon)


###################################################
### code chunk number 445: qt-app-statusbar
###################################################
statusbar <- Qt$QStatusBar()
main_window$setStatusBar(statusbar)


###################################################
### code chunk number 446: qt-app-statusbar-temporary
###################################################
statusbar$showMessage("Load complete", 1000)


###################################################
### code chunk number 447: qt-app-statusbar-perm
###################################################
statusbar$addWidget(Qt$QLabel("Ready"))
statusbar$addPermanentWidget(Qt$QLabel("Version 1.0"))


###################################################
### code chunk number 448: qt-app-dockwidget
###################################################
library(qtutils)
device <- QT()
dock <- Qt$QDockWidget("Device 1")
dock$setWidget(device)


###################################################
### code chunk number 449: qt-app-dock-features
###################################################
dock$features <- Qt$QDockWidget$DockWidgetMovable |
                 Qt$QDockWidget$DockWidgetFloatable


###################################################
### code chunk number 450: qt-app-dock-add
###################################################
main_window$addDockWidget(Qt$Qt$LeftDockWidgetArea, dock)


###################################################
### code chunk number 451: qt-app-dock-tabify
###################################################
device2 <- QT()
dock2 <- Qt$QDockWidget("Device 2", device2)
main_window$tabifyDockWidget(dock, dock2)


###################################################
### code chunk number 452: qt-app-floating
###################################################
dock2$floating <- TRUE


