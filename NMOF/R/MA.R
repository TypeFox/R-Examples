MA <- function(y, order, pad = NULL) {
    order <- as.integer(order)
    if (order < 1L)
        stop("'order' must be a positive integer")
    n <- length(y)
    ma <- cumsum(y) / order
    ma[order:n] <- ma[order:n] - c(0, ma[1L:(n - order)])
    if (!is.null(pad) && order > 1L)
        ma[1L:(order - 1L)] <- pad
    ma
}
