predict.msc.slm <- function (object, newdata, ...) {
    ms = object$ms
    x <- model.matrix(ms, newdata)
    mm <- model.matrix(object, x)
    mm %*% object$slm[[ms$predictLevel]]$lm$coefficients    
}

