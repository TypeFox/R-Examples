d2q <- function(x) {
    if (! is.numeric(x))
        stop("argument must be numeric")
    storage.mode(x) <- "double"
    .Call("d2q", x, PACKAGE = "rcdd")
}
