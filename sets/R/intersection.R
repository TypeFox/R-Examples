set_intersection <-
function(...)
    .set_intersection(list(...))

.set_intersection <-
function(l, matchfun = .exact_match)
{
    len <- length(l)
    if(len < 1L)
        set()
    else if(len < 2L)
        l[[1L]]
    else if(len < 3L) {
        i <- j <- 1
        if (length.set(l[[2L]]) < length.set(l[[1L]])) i <- 2 else j <- 2
        .make_set_from_list(.as.list(l[[i]])[unique(na.omit(matchfun(l[[j]],
                                                                     l[[i]])))])
    } else
        Recall(c(l[1L], list(Recall(l[-1L]))))
}

gset_intersection <-
function(...)
{
    l <- list(...)

    ## compute target support
    support <- do.call(set_intersection, l)

    ## handle empty support
    if (set_is_empty(support))
        return(set())

    ## handle the "ordinary set" case
    if (isTRUE(all(sapply(l, gset_is_set))))
        return(support)

    ## apply connector
    .make_gset_from_list_of_gsets_and_support_and_connector(l, support, .T.)
}

cset_intersection <-
function(...)
{
    l <- list(...)

    ## check matchfun and orderfun
    matchfun <- .check_matchfun(l)
    orderfun <- .check_orderfun(l)

    ## compute target support using correct matchfuns
    support <- cset(.set_intersection(.as.list(l), matchfun),
                    orderfun, matchfun)

    ## handle empty supprt
    if (set_is_empty(support))
        return(set())

    ## handle the "ordinary set" case
    if (isTRUE(all(sapply(l, cset_is_set))))
        return(support)

    ## create gset by applying conorm, and then make cset
    .make_cset_from_gset_and_orderfun_and_matchfun(
         .make_gset_from_list_of_gsets_and_support_and_connector(l,
                                                                 support,
                                                                 .T.,
                                                                 matchfun),
         orderfun,
         matchfun)
}
