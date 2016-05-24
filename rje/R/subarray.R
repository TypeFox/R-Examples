subarray <-
function (x, levels, drop = TRUE) 
{
    if (length(levels) != length(dim(x))) {
        stop("Array and indexlist are not compatible!")
    }
    args <- c(quote(x), levels, list(drop = drop))
    return(do.call("[", args))
}
