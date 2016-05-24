combinations <-
function (p) 
{
    tot = prod(p)
    cp = cumprod(p)
    retval = rep.int(0, tot * length(p))
    dim(retval) = c(tot, length(p))
    for (i in seq_along(p)) {
        retval[, i] <- rep(seq_len(p[i]) - 1, each = cp[i]/p[i])
    }
    retval
}
