`summary.mefa` <-
function(object, nlist = 10, ...)
{
x <- object
out <- list(
    s.rich = rowSums(x$xtab > 0),
    s.abu = rowSums(x$xtab),
    t.occ = colSums(x$xtab > 0),
    t.abu = colSums(x$xtab),
    ntot = sum(x$xtab),
    mfill = sum(x$xtab > 0)/(dim(x)[1]*dim(x)[2]),
    nsamp = dim(x)[1],
    ntaxa = dim(x)[2],
    nsegm = dim(x)[3],
    segment = dimnames(x)$segm,
    call = x$call,
    nested = attr(x, "nested"),
    drop.zero = attr(x, "drop.zero"),
    xtab.fixed = attr(x, "xtab.fixed"))
class(out) <- c("summary.mefa")
attr(out, "nlist") <- nlist
return(out)}

