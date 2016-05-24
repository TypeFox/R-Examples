gset_product <-
function(...)
{
    l <- list(...)

    ## compute target support
    support <- do.call(set_intersection, l)
    if (set_is_empty(support))
        return(set())

    ## handle the "ordinary set" case
    if (isTRUE(all(sapply(l, gset_is_set))))
        return(support)

    ## apply connector
    .make_gset_from_list_of_gsets_and_support_and_connector(l, support, `*`)
}

cset_product <-
function(...)
{
    l <- list(...)

    ## check matchfun and orderfun
    matchfun <- .check_matchfun(l)
    orderfun <- .check_orderfun(l)

    ## compute target support using correct matchfuns
    support <- cset(.set_intersection(.as.list(l), matchfun),
                    orderfun, matchfun)
    if (isTRUE(cset_is_empty(support)))
        return(set())

    ## handle the "ordinary set" case
    if (isTRUE(all(sapply(l, cset_is_set))))
        return(support)

    ## create gset by applying conorm, and then make cset
    .make_cset_from_gset_and_orderfun_and_matchfun(
         .make_gset_from_list_of_gsets_and_support_and_connector(l,
                                                                 support,
                                                                 `*`,
                                                                 matchfun),
         orderfun,
         matchfun)
}
