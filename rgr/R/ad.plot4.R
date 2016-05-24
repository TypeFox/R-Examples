ad.plot4 <-
function (x, xname = deparse(substitute(x)), if.order = FALSE,
    ad.tol = NULL, ifalt = FALSE, ldl = NULL, maxrat = NULL,
    if.text = FALSE, if.cpp = FALSE, ...) 
{
    n <- length(x)
    ndup <- n/2
    x1 <- numeric(ndup)
    x2 <- numeric(ndup)
    if (ifalt) {
        for (i in 1:ndup) {
            j <- 2 * (i - 1) + 1
            x1[i] <- x[j]
            x2[i] <- x[j + 1]
        }
    }
    else {
        for (i in 1:ndup) {
            x1[i] <- x[i]
            x2[i] <- x[ndup + i]
        }
    }
    ad.plot3(x1, x2, xname = xname, if.order = if.order,
            ad.tol = ad.tol, ldl = ldl, maxrat = maxrat, 
            if.text = if.text, if.cpp = if.cpp, ...)
    invisible()
}
