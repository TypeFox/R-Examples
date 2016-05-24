`holdout.Mean.Arima` <-
function(object, test.Vector, ...){
  tmp <- predict(object, n.ahead = length(test.Vector), se.fit = FALSE) - test.Vector
  mean(abs(tmp))
}

