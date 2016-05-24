mse <-
function(yTest, yHat) {
  mean((yTest-yHat)^2)
}
