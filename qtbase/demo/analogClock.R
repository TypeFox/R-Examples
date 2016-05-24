## Draw an analog clock on a widget

## When constructed, starts an update timer
qsetClass("AnalogClock", Qt$QWidget, function(parent = NULL) {
  ## save the timer so that it is not GC'd (and stopped)
  this$timer <- qtimer(1000, update)
  this$timer$start()
  setWindowTitle("Analog Clock")
  resize(200, 200)
})

## Paints the clock
qsetMethod("paintEvent", AnalogClock, function(event) {  
  hourHand <- qpolygon(c(7, -7, 0), c(8, 8, -40))
  minuteHand <- qpolygon(c(7, -7, 0), c(8, 8, -40))

  hourColor <- qcolor(127, 0, 127)
  minuteColor <- qcolor(0, 127, 127, 191)  

  side <- min(c(width, height))

  time <- as.POSIXlt(Sys.time())

  painter <- Qt$QPainter(this)
  painter$setRenderHint(Qt$QPainter$Antialiasing)
  ## coords between -100, 100
  painter$translate(width / 2, height / 2)
  painter$scale(side / 200.0, side / 200.0)

  ## draw the hour hand first
  painter$setPen(Qt$Qt$NoPen)
  painter$setBrush(qbrush(hourColor))
  painter$save()
  painter$rotate(30.0 * ((time$hour + time$min / 60.0)));
  painter$drawConvexPolygon(hourHand)
  painter$restore()

  ## draw hour markers
  painter$setPen(hourColor)

  replicate(12, { 
    painter$drawLine(88, 0, 96, 0)
    painter$rotate(30.0)
  })

  ## draw minute hand
  painter$setPen(Qt$Qt$NoPen)
  painter$setBrush(qbrush(minuteColor))
  painter$save()
  painter$rotate(6.0 * (time$min + time$sec / 60.0))
  painter$drawConvexPolygon(minuteHand)
  painter$restore()

  ## draw minute markers
  painter$setPen(minuteColor)
  sapply(seq(60)-1L, function(i) {
    if ((i %% 5) != 0) # do not draw over the hour markers
      painter$drawLine(92, 0, 96, 0)
    painter$rotate(6.0)
  })

  ## always have to 'end' the painter
  painter$end() 
}, "protected")

clock <- AnalogClock()

print(clock)
