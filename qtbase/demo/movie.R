## Plays animated image formats
qsetClass("MoviePlayer", Qt$QWidget, function(parent = NULL) {
  super(parent)
  
  this$movie <- Qt$QMovie(this)
  movie$setCacheMode(Qt$QMovie$CacheAll)

  this$movieLabel <- Qt$QLabel("No movie loaded")
  movieLabel$setAlignment(Qt$AlignCenter)
  movieLabel$setSizePolicy(Qt$QSizePolicy$Ignored, Qt$QSizePolicy$Ignored)
  movieLabel$setBackgroundRole(Qt$QPalette$Dark)
  movieLabel$setAutoFillBackground(TRUE)

  this$currentMovieDirectory <- "movies"

  createControls()
  createButtons()

  qconnect(movie, "frameChanged", updateFrameSlider)
  qconnect(movie, "stateChanged", updateButtons)
  qconnect(fitCheckBox, "clicked", fitToWindow)
  qconnect(frameSlider, "valueChanged", goToFrame)
  qconnect(speedSpinBox, "valueChanged(int)", movie$setSpeed)
  
  this$mainLayout <- Qt$QVBoxLayout()
  mainLayout$addWidget(movieLabel)
  mainLayout$addLayout(controlsLayout)
  mainLayout$addLayout(buttonsLayout)
  setLayout(mainLayout)

  updateFrameSlider()
  updateButtons()

  setWindowTitle("Movie Player")
  resize(400, 400)
})

qsetMethod("open", MoviePlayer, function() {
  fileName <- Qt$QFileDialog$getOpenFileName(this, "Open a Movie",
                                             currentMovieDirectory)
  if (!is.null(fileName))
    openFile(fileName)
}, "private")

qsetMethod("openFile", MoviePlayer, function(fileName) {
  this$currentMovieDirectory <- Qt$QFileInfo(fileName)$path()

  movie$stop()
  movieLabel$setMovie(movie)
  movie$setFileName(fileName)
  movie$start()

  updateFrameSlider()
  updateButtons()
})

qsetMethod("goToFrame", MoviePlayer, function(frame) {
  movie$jumpToFrame(frame)
}, "private")

qsetMethod("fitToWindow", MoviePlayer, function() {
  movieLabel$setScaledContents(fitCheckBox$isChecked())
}, "private")

qsetMethod("updateFrameSlider", MoviePlayer,
           function(frame = movie$currentFrameNumber())
{
  hasFrames <- frame >= 0

  if (hasFrames) {
    if (movie$frameCount() > 0) {
      frameSlider$setMaximum(movie$frameCount() - 1)
    } else {
      if (frame > frameSlider$maximum())
        frameSlider$setMaximum(frame)
    }
    frameSlider$setValue(frame)
  } else {
    frameSlider$setMaximum(0)
  }
  frameLabel$setEnabled(hasFrames)
  frameSlider$setEnabled(hasFrames)
}, "private")

qsetMethod("updateButtons", MoviePlayer, function(state = movie$state()) {
  playButton$setEnabled(movie$isValid() && movie$frameCount() != 1
                        && state == Qt$QMovie$NotRunning)
  pauseButton$setEnabled(state != Qt$QMovie$NotRunning)
  pauseButton$setChecked(state == Qt$QMovie$Paused)
  stopButton$setEnabled(state != Qt$QMovie$NotRunning)
}, "private")

qsetMethod("createControls", MoviePlayer, function() {
  this$fitCheckBox <- Qt$QCheckBox("Fit to Window")

  this$frameLabel <- Qt$QLabel("Current frame:")

  this$frameSlider <- Qt$QSlider(Qt$Qt$Horizontal)
  frameSlider$setTickPosition(Qt$QSlider$TicksBelow)
  frameSlider$setTickInterval(10)

  this$speedLabel <- Qt$QLabel("Speed:")

  this$speedSpinBox <- Qt$QSpinBox()
  speedSpinBox$setRange(1, 9999)
  speedSpinBox$setValue(100)
  speedSpinBox$setSuffix("%")

  this$controlsLayout <- Qt$QGridLayout()
  controlsLayout$addWidget(fitCheckBox, 0, 0, 1, 2)
  controlsLayout$addWidget(frameLabel, 1, 0)
  controlsLayout$addWidget(frameSlider, 1, 1, 1, 2)
  controlsLayout$addWidget(speedLabel, 2, 0)
  controlsLayout$addWidget(speedSpinBox, 2, 1)
}, "private")

qsetMethod("createButtons", MoviePlayer, function() {
  iconSize <- qsize(36L, 36L)
  
  this$openButton <- Qt$QToolButton()
  openButton$setIcon(style()$standardIcon(Qt$QStyle$SP_DialogOpenButton))
  openButton$setIconSize(iconSize)
  openButton$setToolTip("Open File")
  qconnect(openButton, "clicked", open)

  this$playButton <- Qt$QToolButton()
  playButton$setIcon(style()$standardIcon(Qt$QStyle$SP_MediaPlay))
  playButton$setIconSize(iconSize)
  playButton$setToolTip("Play")
  qconnect(playButton, "clicked", movie$start)
  
  this$pauseButton <- Qt$QToolButton()
  pauseButton$setCheckable(TRUE)
  pauseButton$setIcon(style()$standardIcon(Qt$QStyle$SP_MediaPause))
  pauseButton$setIconSize(iconSize)
  pauseButton$setToolTip("Pause")
  qconnect(pauseButton, "clicked(bool)", movie$setPaused)

  this$stopButton <- Qt$QToolButton()
  stopButton$setIcon(style()$standardIcon(Qt$QStyle$SP_MediaStop))
  stopButton$setIconSize(iconSize)
  stopButton$setToolTip("Stop")
  qconnect(stopButton, "clicked", movie$stop)

  this$quitButton <- Qt$QToolButton()
  quitButton$setIcon(style()$standardIcon(Qt$QStyle$SP_DialogCloseButton))
  quitButton$setIconSize(iconSize)
  quitButton$setToolTip("Quit")
  qconnect(quitButton, "clicked", close)
  
  this$buttonsLayout <- Qt$QHBoxLayout()
  buttonsLayout$addStretch()
  buttonsLayout$addWidget(openButton)
  buttonsLayout$addWidget(playButton)
  buttonsLayout$addWidget(pauseButton)
  buttonsLayout$addWidget(stopButton)
  buttonsLayout$addWidget(quitButton)
  buttonsLayout$addStretch()
}, "private")

MoviePlayer()$show()
