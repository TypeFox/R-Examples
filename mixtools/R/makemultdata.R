makemultdata = function (..., cuts) 
{
    temp = sapply(list(...), length)
    m = max(temp)
    g <- length(cuts) + 1
    cuts <- sort(cuts)
    if (sum(temp != m) > 0) {
        full.data <- function(x, maxm) {
            if (!missing(x)) {
                if (length(x) < maxm) {
                  x <- c(x, NA * rep(1, maxm - length(x)))
                }
            }
            else {
                x <- numeric(0)
            }
            x
        }
        x = sapply(list(...), full.data, maxm = m)
    }
    else {
        if (sapply(list(...), is.matrix) == 1 | sapply(list(...), is.data.frame) == 1) {
            x = t(...)
        }
        else x = cbind(...)
    }
    cutfunc <- function(x, lcut, ucut) {
        x <- na.omit(x)
        sum((x <= ucut) * (x > lcut))
    }
    n <- ncol(x)
    y <- matrix(0, g, n)
    y[1, ] <- apply(x, 2, cutfunc, ucut = cuts[1], lcut = -Inf)
    y[g, ] <- apply(x, 2, cutfunc, ucut = Inf, lcut = cuts[g - 
        1])
    if (g > 2) {
        for (i in 2:(g - 1)) {
            y[i, ] <- apply(x, 2, cutfunc, ucut = cuts[i], lcut = cuts[i - 
                1])
        }
    }
    list(x = t(x), y = t(y))
}