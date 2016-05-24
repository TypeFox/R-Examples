`holdout.Med.lm` <-
function(object, test.Frame, test.Vector, ...){
    tmp <- predict(object, newdata = test.Frame) - test.Vector
    median(abs(tmp))
}

