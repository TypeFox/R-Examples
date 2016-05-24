model.matrix.intReg <- function (object, ...) {
    if (n_match <- match("x", names(object), 0)) 
        object[[n_match]]
    else {
        data <- model.frame(object, xlev = object$xlevels, ...)
        NextMethod("model.matrix", data = data, contrasts = object$contrasts)
    }
}
