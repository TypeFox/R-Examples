image.raw <- function(x, col = c("white", "red", "blue", "green",
    "black", "cyan", "yellow", "magenta"), ...)
{
    stopifnot(is.raw(x))
    foo <- unpackPotts(x)
    bar <- seq(1:nrow(foo))
    baz <- seq(1:ncol(foo))
    qux <- unique(as.vector(foo))
    if (length(qux) > length(col)) stop("not enough colors specified")
    image(bar, baz, foo, axes = FALSE, asp = 1, xlab = "", ylab = "",
       col = col, ...)
}
