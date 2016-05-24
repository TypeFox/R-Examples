`subarray<-` <-
function (x, levels, value) 
{
    if (length(levels) != length(dim(x))) {
        stop("Array and indexlist are not compatible!")
    }
    args <- c(quote(x), levels, quote(value))
    return(do.call("[<-", args))
}
