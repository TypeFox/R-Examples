mape <-
function(yTest, yHat) {
    mean(abs((yTest-yHat)/yTest))*100
  }
