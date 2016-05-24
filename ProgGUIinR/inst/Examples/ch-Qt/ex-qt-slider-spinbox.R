

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

spinbox$value <- slider$value

w$show()
w$raise()

