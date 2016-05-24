confusion.list <-
function (object, truth, ...)
{
    dd <- names(object)
    y <- seq(dd)
    x <- attr(object, "dimension")
    if (!length(x)) 
        x <- seq(dd)
    for (i in y) {
        confi <- confusion(object[, i], truth)
        y[i] <- attr(confi, "error")
    }
    return(x, y)
}

