plot.mritc <- function(x, ...){
    class <- max.col(x$prob)
    x$mask[x$mask==1] <- class
    slices3d(x$mask, ...)
}
