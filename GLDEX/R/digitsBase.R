"digitsBase" <-
function (x, base = 2, ndigits = 1 + floor(log(max(x), base))) 
{
    if (any((x <- as.integer(x)) < 0)) 
        stop("`x' must be non-negative integers")
    r <- matrix(0, nrow = ndigits, ncol = length(x))
    if (ndigits >= 1) 
        for (i in ndigits:1) {
            r[i, ] <- x%%base
            if (i > 1) 
                x <- x%/%base
        }
    class(r) <- "basedInt"
    attr(r, "base") <- base
    r
}

