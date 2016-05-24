`holdout.Mean.lm` <-
function(object, test.Frame, test.Vector, ...){
    tmp <- predict(object, newdata = test.Frame) - test.Vector
    mean(abs(tmp)) 
}

