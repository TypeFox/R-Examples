## Widget for choosing characters. Demonstrates interaction with custom drawing.

qsetClass("CharacterWidget", Qt$QWidget, function(parent = NULL) {
  super(parent)
  this$squareSize <- 24L # size of each character square
  this$columns <- 16L # dimensions of the grid
  this$lastKey <- -1L # no key selected
  this$displayFont <- qfont() # default font
  setMouseTracking(TRUE) # so we receive mouse movement events 
})

## The signal emitted when the selected character changes
qsetSignal("characterSelected(QString character)", CharacterWidget)

## Allow user to change the font family and size

qsetMethod("updateFont", CharacterWidget, function(font) {
  displayFont$setFamily(font$family())
  this$squareSize <- max(c(24L, Qt$QFontMetrics(displayFont)$xHeight() * 3L))
  adjustSize()
  update()
})

qsetMethod("updateSize", CharacterWidget, function(fontSize) {
  displayFont$setPointSize(as.integer(fontSize))
  this$squareSize <- max(c(24L, Qt$QFontMetrics(displayFont)$xHeight() * 3L))
  adjustSize()
  update()
})

qsetMethod("updateStyle", CharacterWidget, function(fontStyle) {
  fontDatabase <- Qt$QFontDatabase()
  oldStrategy <- displayFont$styleStrategy()
  this$displayFont <- fontDatabase$font(displayFont$family(), fontStyle,
                                        displayFont$pointSize())
  displayFont$setStyleStrategy(oldStrategy)
  this$squareSize <- max(c(24L, Qt$QFontMetrics(displayFont)$xHeight() * 3L))
  adjustSize()
  update()
})

qsetMethod("updateFontMerging", CharacterWidget, function(enable) {
  if (enable)
    displayFont$setStyleStrategy(Qt$QFont$PreferDefault)
  else
    displayFont$setStyleStrategy(Qt$QFont$NoFontMerging)
  adjustSize()
  update()
})

## Size is fixed to fit all possible squares
qsetMethod("sizeHint", CharacterWidget, function() {
  qsize(columns*squareSize, (65536L%/%columns)*squareSize)
})

## Paint the characters

qsetMethod("paintEvent", CharacterWidget, function(event) {
  painter <- Qt$QPainter(this)
  redrawRect <- event$rect()  
  painter$fillRect(redrawRect, qbrush("white")) # clear area
  painter$setFont(displayFont)

  ## determine squares to redraw
  beginRow <- redrawRect$top() %/% squareSize
  endRow <- redrawRect$bottom() %/% squareSize
  beginColumn <- redrawRect$left() %/% squareSize
  endColumn <- redrawRect$right() %/% squareSize

  ## draw square borders
  painter$setPen(qpen("gray"))
  sapply(beginRow:endRow, function(row) {
    sapply(beginColumn:endColumn, function(column) {
      painter$drawRect(column*squareSize, row*squareSize,
                       squareSize, squareSize)
    })
  })

  fontMetrics <- Qt$QFontMetrics(displayFont)
  painter$setPen(qpen("black"))

  sapply(beginRow:endRow, function(row) {
    sapply(beginColumn:endColumn, function(column) {
      key <- row*columns + column
      painter$setClipRect(column*squareSize, row*squareSize,
                          squareSize, squareSize)
      if (key == lastKey)
        painter$fillRect(column*squareSize + 1L, row*squareSize + 1L,
                         squareSize, squareSize, qbrush("red"))

      char <- Qt$QChar(key)
      painter$drawText(column*squareSize + (squareSize %/% 2) -
                       fontMetrics$width(char) %/% 2,
                       row*squareSize + 4L + fontMetrics$ascent(),
                       as.character(char))
    })
  })

  painter$end()
}, "protected")

## Respond to mouse press

qsetMethod("mousePressEvent", CharacterWidget, function(event) {
  if (event$button() == Qt$Qt$LeftButton) {
    this$lastKey <- (event$y()%/%squareSize)*columns + event$x()%/%squareSize
    if (Qt$QChar(lastKey)$category() != Qt$QChar$NoCategory)
      characterSelected(as.character(Qt$QChar(lastKey)))
    update()
  }
  else super("mousePressEvent", event)
}, "protected")

## Show tooltip on mouse over

qsetMethod("mouseMoveEvent", CharacterWidget, function(event) {
  widgetPosition <- mapFromGlobal(event$globalPos())
  key <- (widgetPosition$y() %/% squareSize)*columns +
    widgetPosition$x() %/% squareSize;
  text <-
    paste("<p>Character: <span style=\"font-size: 24pt; font-family: ",
          displayFont$family(), "\">", Qt$QChar(key), "</span><p>Value: 0x",
          sprintf("%x", key), sep = "")
  Qt$QToolTip$showText(event$globalPos(), text, this)
}, "protected")

## Main window

qsetClass("MainWindow", Qt$QMainWindow, function() 
          {
            ## Create the widgets
            centralWidget <- Qt$QWidget()
            fontLabel <- Qt$QLabel("Font:")
            this$fontCombo <- Qt$QFontComboBox()
            sizeLabel <- Qt$QLabel("Size:")
            this$sizeCombo <- Qt$QComboBox()
            styleLabel <- Qt$QLabel("Style:")
            this$styleCombo <- Qt$QComboBox()
            fontMergingLabel <- Qt$QLabel("Automatic Font Merging:")
            this$fontMerging <- Qt$QCheckBox()
            fontMerging$setChecked(TRUE)

            this$scrollArea <- Qt$QScrollArea()
            this$characterWidget <- CharacterWidget()
            scrollArea$setWidget(characterWidget)
            
            this$lineEdit <- Qt$QLineEdit()
            clipboardButton <- Qt$QPushButton("&To clipboard")
            this$clipboard <- Qt$QApplication$clipboard()
            
            ## Connect signals
            qconnect(fontCombo, "currentFontChanged", findStyles)
            qconnect(fontCombo, "currentFontChanged", findSizes)
            qconnect(fontCombo, "currentFontChanged",
                     characterWidget$updateFont)
            qconnect(sizeCombo, "currentIndexChanged(QString)",
                     characterWidget$updateSize)
            qconnect(styleCombo, "currentIndexChanged(QString)",
                     characterWidget$updateStyle)
            qconnect(characterWidget, "characterSelected", insertCharacter)
            qconnect(clipboardButton, "clicked", updateClipboard)
            qconnect(fontMerging, "toggled", characterWidget$updateFontMerging)

            ## sync the character widget font with our combo boxes
            characterWidget$updateFont(fontCombo$currentFont)
            findStyles(fontCombo$currentFont)
            findSizes(fontCombo$currentFont)
            
            ## Populate layout
            controlsLayout <- Qt$QHBoxLayout()
            controlsLayout$addWidget(fontLabel)
            controlsLayout$addWidget(fontCombo, 1)
            controlsLayout$addWidget(sizeLabel)
            controlsLayout$addWidget(sizeCombo, 1)
            controlsLayout$addWidget(styleLabel)
            controlsLayout$addWidget(styleCombo, 1)
            controlsLayout$addWidget(fontMergingLabel)
            controlsLayout$addWidget(fontMerging, 1)
            controlsLayout$addStretch(1)

            lineLayout <- Qt$QHBoxLayout()
            lineLayout$addWidget(lineEdit, 1)
            lineLayout$addSpacing(12)
            lineLayout$addWidget(clipboardButton)

            centralLayout <- Qt$QVBoxLayout()
            centralLayout$addLayout(controlsLayout)
            centralLayout$addWidget(scrollArea, 1)
            centralLayout$addSpacing(4)
            centralLayout$addLayout(lineLayout)
            centralWidget$setLayout(centralLayout)

            setCentralWidget(centralWidget)
            setWindowTitle("Character Map")
          })

## Utility functions for populating combo boxes for given font
qsetMethod("findStyles", MainWindow, function(font) {
  currentItem <- styleCombo$currentText
  styleCombo$clear()

  fontDatabase <- Qt$QFontDatabase()
  
  styleCombo$addItems(fontDatabase$styles(font$family()))
  
  styleIndex <- styleCombo$findText(currentItem)

  if (styleIndex == -1)
    styleCombo$setCurrentIndex(0)
  else
    styleCombo$setCurrentIndex(styleIndex)
})
qsetMethod("findSizes", MainWindow, function(font) {
  currentSize <- sizeCombo$currentText
  sizeCombo$blockSignals(TRUE)
  sizeCombo$clear()

  fontDatabase <- Qt$QFontDatabase()
  
  if(fontDatabase$isSmoothlyScalable(font$family(),
                                     fontDatabase$styleString(font)))
    {
      sizeCombo$setEditable(TRUE)
      sizes <- Qt$QFontDatabase$standardSizes() 
    } else {
      sizeCombo$setEditable(FALSE)
      sizes <- fontDatabase$smoothSizes(font$family(),
                                        fontDatabase$styleString(font))
    }

  sizeCombo$addItems(as.character(sizes))
  
  sizeCombo$blockSignals(FALSE)

  sizeIndex <- sizeCombo$findText(currentSize)

  if(sizeIndex == -1)
    sizeCombo$setCurrentIndex(max(c(0, sizeCombo$count / 3)))
  else
    sizeCombo$setCurrentIndex(sizeIndex)
})

## Some signal handlers
qsetMethod("insertCharacter", MainWindow, function(character) {
  lineEdit$insert(character)
})

qsetMethod("updateClipboard", MainWindow, function() {
  clipboard$setText(lineEdit$text, Qt$QClipboard$Clipboard)
  clipboard$setText(lineEdit$text, Qt$QClipboard$Selection)
})

MainWindow()$show()
