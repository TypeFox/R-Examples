qtimer <- function(delay, handler) {
  timer <- Qt$QTimer()
  qconnect(timer, "timeout", handler)
  timer$interval <- delay
  timer
}
