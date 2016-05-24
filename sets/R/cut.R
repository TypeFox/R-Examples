cut.gset <-
function(x, level = 1, type = c("alpha", "nu"), strict = FALSE, ...)
{
    type <- match.arg(type)

    support <- .as.list(.get_support(x))
    memberships <- .get_memberships(x)

    if (type == "alpha") {
        if (gset_is_crisp(x, na.rm = TRUE))
            return(x)

        if (!gset_is_fuzzy_multiset(x)) {
            ind <- if (strict)
                memberships > level
            else
                memberships >= level

            if (length(ind) > 0L)
                .make_gset_from_support_and_memberships(support, as.double(ind))
            else
                set()
        } else {
            M <- sapply(memberships, if (strict)
                        function(i) sum(as.double(unlist(i) > level) * .get_memberships(i))
            else
                        function(i) sum(as.double(unlist(i) >= level) * .get_memberships(i))
                        )
            ind <- is.na(M) | M > 0
            .make_gset_from_support_and_memberships(support[ind], M[ind])
        }
    } else {
        if (gset_is_set_or_fuzzy_set(x, na.rm = TRUE)) {
            if (level > 1)
                set()
            else
                x
        } else if (!gset_is_fuzzy_multiset(x)) {
            ind <- if (strict)
                memberships > level
            else
                memberships >= level

            if (length(ind) > 0L)
                .make_gset_from_support_and_memberships(support, as.double(ind))
            else
                set()

        } else {
            M <- sapply(memberships, function(i) {
                m <- .expand_membership(i)
                if (isTRUE(m[level + strict] == 1))
                    1
                else if (any(is.na(m)))
                    NA
                else if (length(m) < level + strict)
                    0
                else m[level + strict]
            })
            .make_gset_from_support_and_memberships(support, M)
        }
    }
}

cut.cset <-
function(x, level = 1, type = c("alpha", "nu"), strict = FALSE, ...)
    cut.gset(x, level = level, ...)
