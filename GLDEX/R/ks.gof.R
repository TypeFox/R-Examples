ks.gof <-
function (x, y, ..., alternative = c("two.sided", "less", "greater"), 
    exact = NULL) 
{
    if (any(duplicated(x))) {
        dup.ind.x <- duplicated(x)
        x[dup.ind.x] <- jitter(x[dup.ind.x])
    }
    if (!missing(y)) {
        if (any(duplicated(y))) {
            dup.ind.y <- duplicated(y)
            y[dup.ind.y] <- jitter(y[dup.ind.y])
        }
    }
    RVAL <- ks.test(x, y, ..., alternative = c("two.sided", "less", 
        "greater"), exact = NULL)
    return(RVAL)
}
