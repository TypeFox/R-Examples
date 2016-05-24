predict.msr.elnet <- function (object, newdata, ...) {
    ms = object$ms
    x <- model.matrix(ms, newdata)
    mm <- model.matrix(object, x)
    predict(object$slm, mm)
}

