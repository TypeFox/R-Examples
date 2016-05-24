caterpoints <- function(x, parnames, horizontal=TRUE, ...){
    if (!missing(parnames))
        x <- x[parnames]
    if (horizontal)
        points(x, rev(seq(along=x)), ...)
    else
        points(seq(along=x), x, ...)
}
