.sweep0 <-
function (x, MARGIN, STATS, FUN = "-") 
{
    FUN <- match.fun(FUN)
    dims <- dim(x)
    perm <- c(MARGIN, seq_along(dims)[-MARGIN])
    FUN(x, aperm.default(array(STATS, dims[perm]), order(perm)))
}
