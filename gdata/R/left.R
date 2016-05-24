left  <- function(x, n=6L) UseMethod("left")
right <- function(x, n=6L) UseMethod("left")

left.data.frame <- function(x, n=6)
{
    stopifnot(length(n) == 1L)
    n <- if (n < 0L)
        max(ncol(x) + n, 0L)
    else min(n, ncol(x))
    x[, seq_len(n), drop = FALSE]
}
left.matrix <- left.data.frame


right.data.frame <- function (x, n = 6L, ...)
{
    stopifnot(length(n) == 1L)
    ncx <- ncol(x)
    n <- if (n < 0L)
        max(ncx + n, 0L)
    else min(n, ncx)
    x[, seq.int(to = ncx, length.out = n), drop = FALSE]
}

right.matrix <- function (x, n = 6L, addcolnums = TRUE, ...)
{
    stopifnot(length(n) == 1L)
    ncx <- ncol(x)
    n <- if (n < 0L)
        max(ncx + n, 0L)
    else min(n, ncx)
    sel <- seq.int(to = ncx, length.out = n)
    ans <- x[, sel, drop = FALSE]
    if (addcolnums && is.null(colnames(x)))
        colnames(ans) <- paste0("[", sel, ",]")
    ans
}

