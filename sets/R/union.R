c.set <-
set_union <-
function(...)
{
    ret <- do.call(.set,
                   unlist(lapply(list(...), .make_list_elements), recursive = FALSE))
    n <- names(ret)
    if (!is.null(n) && any(duplicated(n[n != ""])))
        names(ret) <- NULL
    ret
}

c.gset <-
gset_union <-
function(...)
{
    l <- list(...)

    ## compute target support
    support <- do.call(set_union, l)

    ## handle the "ordinary set" case
    if (isTRUE(all(sapply(l, gset_is_set))))
        return(support)

    ## apply connector
    .make_gset_from_list_of_gsets_and_support_and_connector(l, support, .S.)
}

c.cset <-
cset_union <-
function(...)
{
    l <- list(...)

    ## check matchfun and orderfun
    matchfun <- .check_matchfun(l)
    orderfun <- .check_orderfun(l)

    ## compute target support using correct matchfuns
    support <- cset(do.call(set_union, l),
                    matchfun = matchfun,
                    orderfun = orderfun)

    ## handle the "ordinary set" case
    if (isTRUE(all(sapply(l, gset_is_set))))
        return(support)

    ## create gset by applying conorm, and then make cset
    .make_cset_from_gset_and_orderfun_and_matchfun(
         .make_gset_from_list_of_gsets_and_support_and_connector(l,
                                                                 support,
                                                                 .S.,
                                                                 matchfun),
         orderfun,
         matchfun)
}

