CloseCommanderRestart <- function() {
  closeCommander(ask.save=getRcmdr("ask.on.exit"))
  Commander()
}

CloseCommanderNoQuestionRestart <- function() {
  closeCommander(ask=FALSE)
  Commander()
}
