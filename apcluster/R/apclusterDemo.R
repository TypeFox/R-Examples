apclusterDemo <- function(l=100, d=2, seed=NA, ...)
{
    if (!is.na(seed)) set.seed(seed)

    if (round(l) != l || l < 2)
        stop("'l' must be an integer at least as large as 2")
    else if (round(d) != d || d < 1)
        stop("'d' must be an integer at least as large as 1")

    x <- matrix(runif(l * d), c(l, d))
    s <- negDistMat(x, r=2)

    # Call function apcluster(), turn on details for later plotting
    apresultObj <- apcluster(s, details=TRUE, ...)

    show(apresultObj)
    plot(apresultObj)

    invisible(list(x, s, apresultObj))
}
