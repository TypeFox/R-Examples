gset_sum <-
function(...)
{
    l <- list(...)

    ## determine connector
    CON <- if (all(sapply(l, gset_is_crisp, na.rm = TRUE)))
        `+`
    else
        function(x, y) pmin(1, x + y)

    ## compute target support
    support <- do.call(set_union, l)

    ## apply connector
    .make_gset_from_list_of_gsets_and_support_and_connector(l, support, CON)

}

gset_difference <-
function(...)
{
    l <- list(...)

    ## connector
    CON <- function(x, y) pmax(0, x - y)

    ## compute target support
    support <- as.set(l[[1]])
    if (set_is_empty(support))
        return(set())

    ## apply connector
    .make_gset_from_list_of_gsets_and_support_and_connector(l, support, CON)
}

cset_sum <-
function(...)
{
    l <- list(...)

    ## check matchfun and orderfun
    matchfun <- .check_matchfun(l)
    orderfun <- .check_orderfun(l)

    ## determine connector
    CON <- if (all(sapply(l, gset_is_crisp, na.rm = TRUE)))
        `+`
    else
        function(x, y) pmin(1, x + y)

    ## compute target support
    support <- cset(do.call(set_union, l),
                    matchfun = matchfun,
                    orderfun = orderfun)

    ## apply connector
    .make_cset_from_gset_and_orderfun_and_matchfun(
        .make_gset_from_list_of_gsets_and_support_and_connector(l,
                                                                support,
                                                                CON,
                                                                matchfun),
        orderfun,
        matchfun)

}

cset_difference <-
function(...)
{
    l <- list(...)

    ## check matchfun and orderfun
    matchfun <- .check_matchfun(l)
    orderfun <- .check_orderfun(l)

    ## connector
    CON <- function(x, y) pmax(0, x - y)

    ## compute target support
    support <- cset(as.set(l[[1]]),
                    matchfun = matchfun,
                    orderfun = orderfun)
    if (isTRUE(cset_is_empty(support)))
        return(support)


    ## apply connector
    .make_cset_from_gset_and_orderfun_and_matchfun(
        .make_gset_from_list_of_gsets_and_support_and_connector(l,
                                                                support,
                                                                CON,
                                                                matchfun),
        orderfun,
        matchfun)
}

