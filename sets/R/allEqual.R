all.equal.set <-
function(target, current, ...)
{
    if (gset_is_equal(target, current, na.rm = TRUE))
        return(TRUE)

    msg <- NULL
    if (data.class(target) != data.class(current)) {
        msg <- paste("target is ", data.class(target), ", current is ",
                     data.class(current), sep = "")
    } else {
        if (isTRUE((lt <- length(target)) != (lc <- length(current))))
            msg <- paste("cardinality of target is ", lt, ", of current is ",
                         lc, sep = "")
        if (isTRUE((lc <- length(gset_difference(current, target))) > 0))
            msg <- c(msg, paste(lc, "elements of current are not in target"))
        if (isTRUE((lt <- length(gset_difference(target, current))) > 0))
            msg <- c(msg, paste(lt, "elements of target are not in current"))
    }
    msg
}

all.equal.gset <-
function(target, current, ...)
{
    if (gset_is_equal(target, current, na.rm = TRUE))
        return(TRUE)
    if (isTRUE(all.equal(as.list(target), as.list(current)))
        && isTRUE(all.equal(gset_memberships(target), gset_memberships(current))))
        return(TRUE)

    if (data.class(target) != data.class(current))
        paste("target is ", data.class(target), ", current is ",
              data.class(current), sep = "")
    else if (gset_is_fuzzy_set(target, na.rm = TRUE) &&
             gset_is_multiset(current, na.rm = TRUE))
        "target is a fuzzy set, current is a multiset"
    else if (gset_is_fuzzy_set(current, na.rm = TRUE) &&
             gset_is_multiset(target, na.rm = TRUE))
        "current is a fuzzy set, target is a multiset"
    else if (!isTRUE(all.equal(lt <- gset_cardinality(target, na.rm = TRUE),
             lc <- gset_cardinality(current, na.rm = TRUE))))
        paste("cardinality of target is ", lt, ", of current is ",
              lc, sep = "")
    else
        "target and current have different memberships"
}

all.equal.cset <-
function(target, current, ...)
{
    if (cset_is_equal(target, current, na.rm = TRUE))
        return(TRUE)

    msg <- all.equal.gset(target, current, ...)
    if (isTRUE(msg)) msg <- ""
    if (!identical(.orderfun(target), .orderfun(current)))
        msg <- c(msg, "order functions of target and current differ")
    if (!identical(.matchfun(target), .matchfun(current)))
        msg <- c(msg, "match functions of target and current differ")
    msg
}
