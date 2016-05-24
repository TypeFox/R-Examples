predict.msr <- function (object, newdata, ...) {
    ms = object$ms
    x <- model.matrix(ms, newdata)
    mm <- model.matrix(object, x)
    mm %*% object$slm$lm$coefficients    
}

