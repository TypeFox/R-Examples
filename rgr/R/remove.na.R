remove.na <-
function (xx, iftell = TRUE) 
{
    if (is.vector(xx)) {
        nx <- length(xx)
        x <- na.omit(xx)
        n <- length(x)
        m <- 1
        intype <- "vector"
    }
    else {
        nx <- length(xx[, 1])
        x <- na.omit(xx)
        n <- length(x[, 1])
        m <- length(x[1, ])
        intype <- "matrix"
    }
    nna <- nx - n
    if (iftell & nna > 0) {
        if (intype == "matrix") 
            cat(paste("\n ", nna, "row(s) with NA(s) removed from matrix\n"))
        else {
            cat(paste("\n ", nna, "NA(s) removed from vector\n"))
        }
    }
    invisible(list(x = x, n = n, m = m, nna = nna))
}
