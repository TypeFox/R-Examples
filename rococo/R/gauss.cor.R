gauss.cor <- function (x, y)
{
    if (!is.vector(x) || !is.vector(y) || length(x) != length(y))
        stop("'x' and 'y' must be vectors with the same length")

    cor(qnorm(rank(x) / (length(x) + 1)), qnorm(rank(y) / (length(y) + 1)),
        method="pearson")
}
