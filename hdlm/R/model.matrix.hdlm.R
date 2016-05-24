model.matrix.hdlm <-
function(object, ...)
{
    if(n_match <- match("x", names(object), 0L)) object[[n_match]]
    else {
        data <- model.frame(object, xlev = object$xlevels, ...)
        NextMethod("model.matrix", data = data,
                   contrasts.arg = object$contrasts)
    }
}

