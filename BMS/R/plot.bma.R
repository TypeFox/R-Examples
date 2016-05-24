plot.bma <-
function (x, ...) 
{
    if (!is.bma(x)) 
        stop("Need to provide object of class 'bma'!")
    if (x$arguments$nmodel < 3) {
        try(plotModelsize(x, ...), silent = TRUE)
    }
    else {
        layout(matrix(1:2, 2, 1))
        try(plotModelsize(x, ...), silent = TRUE)
        try(plotConv(x, ...), silent = TRUE)
        layout(1)
    }
}
