rep.set <-
function(x, ...)
    x

rep.gset <-
function(x, ...)
{
    l <- list(...)
    each <- if (length(l) == 0L)
        1
    else
        l[[1]]
    m <- .get_memberships(x)
    memberships <- if (gset_is_crisp(x, na.rm = TRUE))
        m * each
    else if (gset_is_fuzzy_set(x, na.rm = TRUE))
        lapply(m, rep, each = each)
    else
        lapply(m, function(i) gset(.get_support(i), .get_memberships(i) * each))

    .make_gset_from_support_and_memberships(.get_support(x),
                                           .canonicalize_memberships(memberships))
}

rep.cset <-
function(x, ...)
    rep.gset(x, ...)
