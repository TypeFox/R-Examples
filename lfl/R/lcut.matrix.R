lcut3.matrix <- function(x, ...) {
    if (!is.matrix(x)) {
        stop("'x' must be a matrix")
    }
    result <- lcut3(as.data.frame(x), ...)
    return(result)
}


lcut5.matrix <- function(x, ...) {
    if (!is.matrix(x)) {
        stop("'x' must be a matrix")
    }
    result <- lcut5(as.data.frame(x), ...)
    return(result)
}
