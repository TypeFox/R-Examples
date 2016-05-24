set_cardinality <-
length.gset <-
length.cset <-
length.set <-
function(x)
    length(.as.list(x))

gset_cardinality <-
function(x, type = c("absolute", "relative"), na.rm = FALSE)
{
    type <- match.arg(type)

    ret <- if (isTRUE(gset_is_set(x)))
        length.set(x)
    else if (gset_is_fuzzy_multiset(x))
        sum(sapply(.get_memberships(x),
                   function(i) sum(as.double(i) * .get_memberships(i),
                                   na.rm = na.rm)),
            na.rm = na.rm)
    else
        sum(.get_memberships(x), na.rm = na.rm)

    if (type == "absolute")
        ret
    else
        ret / length.set(x)
}

cset_cardinality <-
function(x, type = c("absolute", "relative"), na.rm = FALSE)
    gset_cardinality(x, type, na.rm)

