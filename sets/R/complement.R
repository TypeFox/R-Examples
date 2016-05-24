set_complement <-
function(x, y)
    .set_complement_using_matchfun(x, y)

.set_complement_using_matchfun <-
function(x, y, matchfun = .exact_match)
{
    if (missing(y))
        return(set())
    y <- .as.list(y)
    ind <- unique(na.omit(matchfun(x, y)))
    .make_set_from_list(if(length(ind)) y[-ind] else y)
}

gset_complement <-
function(x, y = NULL)
{
    Cx <- .gset_complement(x)
    if (is.null(y)) Cx else gset_intersection(Cx, y)
}

cset_complement <-
function(x, y = NULL)
{
    if (!is.null(y)) {
        matchfun <- .check_matchfun(list(x, y))
        orderfun <- .check_orderfun(list(x, y))
    } else {
        matchfun <- .matchfun(x)
        orderfun <- .orderfun(x)
    }

    Cx <- .make_cset_from_gset_and_orderfun_and_matchfun(
              .gset_complement(x,
                               universe = cset_universe(x),
                               bound = cset_bound(x)),
               orderfun,
               matchfun
    )
    if (is.null(y)) Cx else cset_intersection(Cx, y)
}

.gset_complement <-
function(x, universe = gset_universe(x), bound = gset_bound(x))
{
    if (set_is_empty(universe) || bound < 1L)
        return(set())

    ## efficiency hack: by default, the complement gets the universe
    ## of the original set, *except* if there is a default universe
    ## and the universe attribute of the original set is missing.
    tuniverse <-
        if (!is.null(sets_options("universe")) && is.null(attr(x, "universe")))
            NULL
        else
            universe

    if (isTRUE(gset_is_set(x)) && bound == 1L)
        gset(set_complement(x, universe), universe = tuniverse, bound = 1L)
    else if (gset_is_crisp(x, na.rm = TRUE)) {
        M <- gset(universe, rep(bound, length(universe)))
        .structure(.set_bound(.set_universe(gset_difference(M, x), tuniverse), bound),
                  class = c("gset", "cset"))
    } else if (gset_is_fuzzy_set(x, na.rm = TRUE) && bound == 1L)
        .make_gset_from_support_and_memberships(
              universe,
              .N.(.memberships_for_support(x, universe)),
              universe = tuniverse,
              bound = bound
        )
    else {
        connector <- function(x, y) .N.(y)
        M <- gset(universe, rep(bound, length(universe)))
        memberships <-
            .apply_connector_to_list_of_gsets_and_support(list(M, x),
                                                          universe,
                                                          connector)
        .make_gset_from_support_and_memberships(
              universe,
              memberships,
              universe = tuniverse,
              bound = bound
        )
    }
}
