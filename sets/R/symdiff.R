set_symdiff <-
function(...)
{
    len <- length(l <- list(...))
    SD <- function(x, y) set_union(set_complement(x, y), set_complement(y, x))

    if (len < 1L)
        gset()
    else if (len < 2L)
        l[[1L]]
    else
        Reduce(SD, l)
}

gset_symdiff <-
function(...)
{
    len <- length(l <- list(...))
    SD <- function(x, y) gset_union(gset_complement(x, y), gset_complement(y, x))
    if (len < 1L)
        gset()
    else if (len < 2L)
        l[[1L]]
    else
        Reduce(SD, l)
}

cset_symdiff <-
function(...)
{
    len <- length(l <- list(...))
    SD <- function(x, y) cset_union(cset_complement(x, y), cset_complement(y, x))
    if (len < 1L)
        set()
    else if (len < 2L)
        l[[1L]]
    else
        Reduce(SD, l)
}

"%D%" <-
function(x, y)
{
    if (is.set(x) && is.set(y))
        set_symdiff(x, y)
    else if (is.gset(x) && is.gset(y))
        gset_symdiff(x, y)
    else if (is.cset(x) && is.cset(y))
        cset_symdiff(x, y)
    else
        stop("Operator only defined for (g,c)sets.")
}
