mae <-
function(yTest, yHat) {
  mean(abs(yTest-yHat))
}
