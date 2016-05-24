digitsBase <- function (x, base = 2, ndigits = 1 + floor(log(max(x), base)))
{
     ## from package sfsmisc
     ## as its import causes problems
     ## internal
     ## class omitted
    if (any(x < 0))
        stop("'x' must be non-negative integers")
    if (any(x != trunc(x)))
        stop("'x' must be integer-valued")
    r <- matrix(0, nrow = ndigits, ncol = length(x))
    if (ndigits >= 1)
        for (i in ndigits:1) {
            r[i, ] <- x%%base
            if (i > 1)
                x <- x%/%base
        }
    attr(r, "base") <- base
    r
}
