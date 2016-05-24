mintersect <- function(..., sorted=FALSE) {
    x <- Reduce(intersect, list(...))
    if (sorted) {
        sort(x)
    } else {
        x
    }
}

munion <- function(..., sorted=FALSE) {
    x <- Reduce(union, list(...))
    if (sorted) {
        sort(x)
    } else {
        x
    }
}
