`holdout.SS.default` <-
function(object, test.Frame, test.Vector, ...){
  tmp <- predict(object, n.ahead = length(test.Vector), se.fit = FALSE, 
  newdata = test.Frame) - test.Vector
  crossprod(tmp, tmp)
}

