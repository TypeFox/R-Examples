require(qtbase)

## This shows basic drag and drop framework


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


### Show this off
w <- Qt$QWidget()
lyt <- Qt$QHBoxLayout()
w$setLayout(lyt)
lyt$addWidget(DragLabel("drag me"))
lyt$addWidget(DropLabel("Drop here"))
w$show()
w$raise()
