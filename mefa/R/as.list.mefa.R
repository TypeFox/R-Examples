`as.list.mefa` <-
function(x, ...)
{
    if (is.null(x$segm) || length(x$segm) == 1) {
        return(list(x$xtab))
    } else {
        return(x$segm)
    }
}

