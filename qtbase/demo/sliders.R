## Three types of sliders: slider, scrollbar, dial
## Synchronized with spin buttons

qsetClass("SlidersGroup", Qt$QGroupBox,
          function(orientation, title, parent = NULL)
{
  super(title, parent)

  ## add the three types of slider widgets
  
  this$slider <- Qt$QSlider(orientation)
  slider$setFocusPolicy(Qt$Qt$StrongFocus)
  slider$setTickPosition(Qt$QSlider$TicksBothSides)
  slider$setTickInterval(10)
  slider$setSingleStep(1)

  this$scrollBar <- Qt$QScrollBar(orientation)
  scrollBar$setFocusPolicy(Qt$Qt$StrongFocus)

  this$dial <- Qt$QDial()
  dial$setFocusPolicy(Qt$StrongFocus)

  ## chain their state changes together, eventually emitting our signal
  qconnect(slider, "valueChanged", scrollBar$setValue)
  qconnect(scrollBar, "valueChanged", dial$setValue)
  qconnect(dial, "valueChanged", slider$setValue)
  qconnect(dial, "valueChanged", valueChanged)
  
  if (orientation == Qt$Qt$Horizontal)
    direction <- Qt$QBoxLayout$TopToBottom
  else
    direction <- Qt$QBoxLayout$LeftToRight

  slidersLayout <- Qt$QBoxLayout(direction)
  slidersLayout$addWidget(slider)
  slidersLayout$addWidget(scrollBar)
  slidersLayout$addWidget(dial)
  setLayout(slidersLayout)
})

## Emitted when one of the sliders changes

qsetSignal("valueChanged(int value)", SlidersGroup)

## Change parameters of sliders in batch

qsetMethod("setValue", SlidersGroup, function(value) {
  slider$setValue(value) # initiate the cascade
})

qsetMethod("setMinimum", SlidersGroup, function(value) {
  slider$setMinimum(value)
  scrollBar$setMinimum(value)
  dial$setMinimum(value)
})

qsetMethod("setMaximum", SlidersGroup, function(value) {
  slider$setMaximum(value)
  scrollBar$setMaximum(value)
  dial$setMaximum(value)
})

qsetMethod("invertAppearance", SlidersGroup, function(invert) {
  slider$setInvertedAppearance(invert)
  scrollBar$setInvertedAppearance(invert)
  dial$setInvertedAppearance(invert)
})

qsetMethod("invertKeyBindings", SlidersGroup, function(invert) {
  slider$setInvertedControls(invert)
  scrollBar$setInvertedControls(invert)
  dial$setInvertedControls(invert)
})

qsetClass("Window", Qt$QWidget, function() {
  this$horizontalSliders <- SlidersGroup(Qt$Qt$Horizontal, "Horizontal")
  this$verticalSliders <- SlidersGroup(Qt$Qt$Vertical, "Vertical")

  this$stackedWidget <- Qt$QStackedWidget()
  stackedWidget$addWidget(horizontalSliders)
  stackedWidget$addWidget(verticalSliders)

  createControls("Controls")

  ## chain together the two slider groups and the spin buttons
  qconnect(horizontalSliders, "valueChanged", verticalSliders$setValue)
  qconnect(verticalSliders, "valueChanged", valueSpinBox$setValue)
  qconnect(valueSpinBox, "valueChanged", horizontalSliders$setValue)
  
  layout <- Qt$QHBoxLayout()
  layout$addWidget(controlsGroup)
  layout$addWidget(stackedWidget)
  setLayout(layout)

  minimumSpinBox$setValue(0)
  maximumSpinBox$setValue(20)
  valueSpinBox$setValue(5)

  setWindowTitle("Sliders")
})

qsetMethod("createControls", Window, function(title) {
  ## Construct the controls
  
  this$controlsGroup <- Qt$QGroupBox(title)

  this$minimumLabel <- Qt$QLabel("Minimum value:")
  this$maximumLabel <- Qt$QLabel("Maximum value:")
  this$valueLabel <- Qt$QLabel("Current value:")

  this$invertedAppearance <- Qt$QCheckBox("Inverted appearance")
  this$invertedKeyBindings <- Qt$QCheckBox("Inverted key bindings")

  this$minimumSpinBox <- Qt$QSpinBox()
  minimumSpinBox$setRange(-100, 100)
  minimumSpinBox$setSingleStep(1)

  this$maximumSpinBox <- Qt$QSpinBox()
  maximumSpinBox$setRange(-100, 100)
  maximumSpinBox$setSingleStep(1)

  this$valueSpinBox <- Qt$QSpinBox()
  valueSpinBox$setRange(-100, 100)
  valueSpinBox$setSingleStep(1)

  this$orientationCombo <- Qt$QComboBox()
  orientationCombo$addItem("Horizontal slider-like widgets")
  orientationCombo$addItem("Vertical slider-like widgets")

  ## Update sliders in response to user input
  
  qconnect(orientationCombo, "activated(int)", stackedWidget$setCurrentIndex)
  qconnect(minimumSpinBox, "valueChanged(int)", horizontalSliders$setMinimum)
  qconnect(minimumSpinBox, "valueChanged(int)", verticalSliders$setMinimum)
  qconnect(maximumSpinBox, "valueChanged(int)", horizontalSliders$setMaximum)
  qconnect(maximumSpinBox, "valueChanged(int)", verticalSliders$setMaximum)
  qconnect(invertedAppearance, "toggled", horizontalSliders$invertAppearance)
  qconnect(invertedAppearance, "toggled", verticalSliders$invertAppearance)
  qconnect(invertedKeyBindings, "toggled", verticalSliders$invertKeyBindings)
  qconnect(invertedKeyBindings, "toggled", horizontalSliders$invertKeyBindings)

  controlsLayout <- Qt$QGridLayout()
  controlsLayout$addWidget(minimumLabel, 0, 0)
  controlsLayout$addWidget(maximumLabel, 1, 0)
  controlsLayout$addWidget(valueLabel, 2, 0)
  controlsLayout$addWidget(minimumSpinBox, 0, 1)
  controlsLayout$addWidget(maximumSpinBox, 1, 1)
  controlsLayout$addWidget(valueSpinBox, 2, 1)
  controlsLayout$addWidget(invertedAppearance, 0, 2)
  controlsLayout$addWidget(invertedKeyBindings, 1, 2)
  controlsLayout$addWidget(orientationCombo, 3, 0, 1, 3)
  controlsGroup$setLayout(controlsLayout)
}, "private")

Window()$show()
