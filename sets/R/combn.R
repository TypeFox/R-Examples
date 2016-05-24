set_combn <-
function(x, m)
{
    if (m == 0)
        set()
    else
        do.call(set, apply(combn(.as.list(x), m), 2L, as.set))
}

gset_combn <-
function(x, m)
{
    if (m == 0L)
        gset()
    else {
        support <- apply(combn(.as.list(x), m), 2L, as.set)
        memberships <- unlist(apply(combn(.get_memberships(x), m), 2L, list),
                              recursive = FALSE)
        do.call(set, Map(gset, support, memberships))
    }
}

cset_combn <-
function(x, m)
{
    ## get matchfun and orderfun
    matchfun <- .matchfun(x)
    orderfun <- .orderfun(x)
    if (is.integer(orderfun))
        orderfun <- NULL

    if (m == 0L)
        gset()
    else {
        support <- apply(combn(.as.list(x), m), 2L, as.set)
        memberships <- unlist(apply(combn(.get_memberships(x), m), 2L, list),
                              recursive = FALSE)
        do.call(set,
                lapply(Map(gset, support, memberships),
                       cset, orderfun, matchfun)
                       )
    }
}

