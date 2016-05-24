predict.msc.slm.elnet <- function (object, newdata, ...) {
    ms = object$ms
    x <- model.matrix(ms, newdata) 
    mm <- model.matrix(object, x)
    predict(object$slm[[ms$predictLevel]]$elnet, as.matrix(mm))   
}

