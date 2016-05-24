.meanfun <-
list(arithmetic = function(x, y) (x + y) / 2,
     geometric = function(x, y) sqrt(x * y),
     harmonic = function(x, y) 2 / (1 / x + 1 / y)
     )

gset_mean <-
function(x, y, type = c("arithmetic", "geometric", "harmonic"))
{
    l <- list(x, y)
    FUN <- .meanfun[[match.arg(type)]]

    ## compute target support
    support <- do.call(set_union, l)

    ## apply connector
    .make_gset_from_list_of_gsets_and_support_and_connector(l, support, FUN,
                                                            enforce_general_case = TRUE)
}

cset_mean <-
function(x, y, type = c("arithmetic", "geometric", "harmonic"))
{
    l <- list(x, y)
    FUN <- .meanfun[[match.arg(type)]]

    ## check matchfun and orderfun
    matchfun <- .check_matchfun(l)
    orderfun <- .check_orderfun(l)

    ## compute target support using correct matchfuns
    support <- cset(do.call(set_union, l),
                    matchfun = matchfun,
                    orderfun = orderfun)

    ## create gset by applying conorm, and then make cset
    .make_cset_from_gset_and_orderfun_and_matchfun(
         .make_gset_from_list_of_gsets_and_support_and_connector(l,
                                                                 support,
                                                                 FUN,
                                                                 matchfun,
                                                                 TRUE),
         orderfun,
         matchfun)
}

