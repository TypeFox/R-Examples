`as.data.frame.stcs` <-
function(x, ...)
{
    attr(x, "call") <- NULL
    attr(x, "expand") <- NULL
    attr(x, "zero.count") <- NULL
    attr(x, "zero.pseudo") <- NULL
    class(x) <- "data.frame"
    x
}

