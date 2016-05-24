## Image viewer with zoom/pan and print function

qsetClass("ImageViewer", Qt$QMainWindow, function() {
  ## Qt displays images with QLabel
  this$imageLabel <- Qt$QLabel()
  imageLabel$setBackgroundRole(Qt$QPalette$Base)
  imageLabel$setSizePolicy(Qt$QSizePolicy$Ignored, Qt$QSizePolicy$Ignored)
  imageLabel$setScaledContents(TRUE)

  this$scrollArea <- Qt$QScrollArea()
  scrollArea$setBackgroundRole(Qt$QPalette$Dark)
  scrollArea$setWidget(imageLabel)
  setCentralWidget(scrollArea)

  createActions()
  createMenus()

  this$printer <- Qt$QPrinter()
  
  setWindowTitle("Image Viewer")
  resize(500, 400)
})

qsetMethod("open", ImageViewer, function() {
  fileName <- Qt$QFileDialog$getOpenFileName(this, "Open File",
                                             Qt$QDir$currentPath())
  if (nzchar(fileName)) {
    image <- Qt$QImage(fileName)
    if (image$isNull()) {
      Qt$QMessageBox$information(this, "Image Viewer",
                                 sprintf("Cannot load %s.", fileName))
      return()
    }
    imageLabel$setPixmap(Qt$QPixmap$fromImage(image))
    this$scaleFactor <- 1.0

    printAct$setEnabled(TRUE)
    fitToWindowAct$setEnabled(TRUE)
    updateActions()

    if (!fitToWindowAct$isChecked())
      imageLabel$adjustSize()
  }
}, "private")

qsetMethod("print", ImageViewer, function() {
  dialog <- Qt$QPrintDialog(printer, this)
  if (dialog$exec()) {
    painter <- Qt$QPainter(printer)
    rect <- painter$viewport()
    size <- imageLabel$pixmap$size()
    size$scale(rect$size(), Qt$Qt$KeepAspectRatio)
    painter$setViewport(rect$x(), rect$y(), size$width(), size$height())
    painter$setWindow(imageLabel$pixmap$rect())
    painter$drawPixmap(0, 0, imageLabel$pixmap)
    painter$end()
  }
}, "private")

qsetMethod("zoomIn", ImageViewer, function() {
  scaleImage(1.25)
}, "private")

qsetMethod("zoomOut", ImageViewer, function() {
  scaleImage(0.8)
}, "private")

qsetMethod("normalSize", ImageViewer, function() {
  imageLabel$adjustSize()
  this$scaleFactor <- 1.0
}, "private")

qsetMethod("fitToWindow", ImageViewer, function() {
  fitToWindow <- fitToWindowAct$isChecked()
  scrollArea$setWidgetResizable(fitToWindow)
  if (!fitToWindow) {
    normalSize()
  }
  updateActions()
}, "private")

qsetMethod("about", ImageViewer, function() {
  Qt$QMessageBox$about(this, "About Image Viewer",
                       paste("<p>The <b>Image Viewer</b> example shows how",
                             "to combine QLabel and QScrollArea to display",
                             "an image. QLabel is typically used for",
                             "displaying a text, but it can also display an",
                             "image. QScrollArea provides a scrolling view",
                             "around another widget. If the child widget",
                             "exceeds the size of the frame, QScrollArea",
                             "automatically provides scroll bars. </p><p>The",
                             "example demonstrates how QLabel's ability to",
                             "scale its contents (QLabel::scaledContents), and",
                             "QScrollArea's ability to automatically resize",
                             "its contents (QScrollArea::widgetResizable),",
                             "can be used to implement zooming and scaling",
                             "features. </p><p>In addition the example shows",
                             "how to use QPainter to print an image.</p>"))
}, "private")

qsetMethod("createActions", ImageViewer, function() {
  this$openAct <- Qt$QAction("&Open...", this)
  openAct$setShortcut(Qt$QKeySequence("Ctrl+O"))
  qconnect(openAct, "triggered", open)
  
  this$printAct <- Qt$QAction("&Print...", this)
  printAct$setShortcut(Qt$QKeySequence("Ctrl+P"))
  printAct$setEnabled(FALSE)
  qconnect(printAct, "triggered", print)

  this$exitAct <- Qt$QAction("E&xit", this)
  exitAct$setShortcut(Qt$QKeySequence("Ctrl+Q"))
  qconnect(exitAct, "triggered", close)

  this$zoomInAct <- Qt$QAction("Zoom &In (25%)", this)
  zoomInAct$setShortcut(Qt$QKeySequence("Ctrl++"))
  zoomInAct$setEnabled(FALSE)
  qconnect(zoomInAct, "triggered", zoomIn)
           
  this$zoomOutAct <- Qt$QAction("Zoom &Out (25%)", this)
  zoomOutAct$setShortcut(Qt$QKeySequence("Ctrl+-"))
  zoomOutAct$setEnabled(FALSE)
  qconnect(zoomOutAct, "triggered", zoomOut)
  
  this$normalSizeAct <- Qt$QAction("&Normal Size", this)
  normalSizeAct$setShortcut(Qt$QKeySequence("Ctrl+S"))
  normalSizeAct$setEnabled(FALSE)
  qconnect(normalSizeAct, "triggered", normalSize)
  
  this$fitToWindowAct <- Qt$QAction("&Fit to Window", this)
  fitToWindowAct$setEnabled(FALSE)
  fitToWindowAct$setCheckable(TRUE)
  fitToWindowAct$setShortcut(Qt$QKeySequence("Ctrl+F"))
  qconnect(fitToWindowAct, "triggered", fitToWindow)
  
  this$aboutAct <- Qt$QAction("&About", this)
  qconnect(aboutAct, "triggered", about)
  
  this$aboutQtAct <- Qt$QAction("About &Qt", this)
  qconnect(aboutQtAct, "triggered", Qt$QApplication$aboutQt)
}, "private")

qsetMethod("createMenus", ImageViewer, function() {
  this$fileMenu <- Qt$QMenu("&File", this)
  fileMenu$addAction(openAct)
  fileMenu$addAction(printAct)
  fileMenu$addSeparator()
  fileMenu$addAction(exitAct)

  this$viewMenu <- Qt$QMenu("&View", this)
  viewMenu$addAction(zoomInAct)
  viewMenu$addAction(zoomOutAct)
  viewMenu$addAction(normalSizeAct)
  viewMenu$addSeparator()
  viewMenu$addAction(fitToWindowAct)

  this$helpMenu <- Qt$QMenu("&Help", this)
  helpMenu$addAction(aboutAct)
  helpMenu$addAction(aboutQtAct)

  menuBar()$addMenu(fileMenu)
  menuBar()$addMenu(viewMenu)
  menuBar()$addMenu(helpMenu)
}, "private")

qsetMethod("updateActions", ImageViewer, function() {
  zoomInAct$setEnabled(!fitToWindowAct$isChecked())
  zoomOutAct$setEnabled(!fitToWindowAct$isChecked())
  normalSizeAct$setEnabled(!fitToWindowAct$isChecked())
}, "private")

qsetMethod("scaleImage", ImageViewer, function(factor) {
  this$scaleFactor <- scaleFactor * factor
  imageLabel$resize(scaleFactor * imageLabel$pixmap$size())

  adjustScrollBar(scrollArea$horizontalScrollBar(), factor)
  adjustScrollBar(scrollArea$verticalScrollBar(), factor)

  zoomInAct$setEnabled(scaleFactor < 3.0)
  zoomOutAct$setEnabled(scaleFactor > 0.333)
}, "private")

qsetMethod("adjustScrollBar", ImageViewer, function(scrollBar, factor) {
  scrollBar$setValue(factor * scrollBar$value +
                     ((factor - 1) * scrollBar$pageStep/2))
}, "private")

ImageViewer()$show()
