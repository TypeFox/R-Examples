### Membership transformations

gset_transform_memberships <-
function(x, FUN, ...)
{
    F <- function(x) pmax(0, pmin(1, FUN(x, ...)))
    m <- if (cset_is_set_or_multiset(x, na.rm = TRUE))
        lapply(.get_memberships(x), function(i) gset(F(1), i))
    else if (!cset_is_fuzzy_multiset(x))
        F(.get_memberships(x))
    else
        lapply(.get_memberships(x), F)
    .set_memberships(x, m)
}

cset_transform_memberships <-
function(x, FUN, ...)
    gset_transform_memberships(x, FUN, ...)

cset_concentrate <-
gset_concentrate <-
function(x)
    gset_transform_memberships(x, function(i) i * i)

cset_dilate <-
gset_dilate <-
function(x)
    gset_transform_memberships(x, sqrt)

cset_normalize <-
gset_normalize <-
function(x, height = 1)
{
    if (height < 0 || height > 1)
        stop("Height must be in the unit interval.")
    gset_transform_memberships(x, function(i) height * i /
                               max(.get_memberships(x), na.rm = TRUE))
}

## internal functions to handle memberships

.expand_membership <-
function(x, decreasing = TRUE, len = NA, rep = TRUE)
{
    if (is.null(x))
        x <- 0
    if (!is.atomic(x))
        x <- .as.list(x)
    M <- .get_memberships(x)

    ## create vector from gset
    ret <- if (rep) {
        if (is.list(x))
            rep.int(as.numeric(x), M)
        else if (length(x) == 1L && !is.na(x) && x > 1)
            rep.int(1, x)
        else
            x
    } else {
        as.numeric(x)
    }

    ## optionally, fill up with 0s
    if (!is.na(len)) {
        ret <- c(ret, rep.int(0, len - length(ret)))
        if (!rep)
            M <- c(M, rep.int(0, len - length(M)))
    }

    ## optionally, sort values
    if (!is.null(decreasing)) {
        O <- order(ret, decreasing = decreasing)
        ret <- ret[O]
        if (!rep) M <- M[O]
    }

    if (rep)
        ret
    else
        .structure(ret, memberships = M)
}

.memberships_for_support <-
function(x, support, matchfun = .exact_match)
{
    m <- matchfun(support, .get_support(x))
    tmp <- if (all(is.na(m)))
        m
    else
        .get_memberships(x)[m]
    tmp[is.na(m)] <- 0
    tmp
}

.canonicalize_memberships <-
function(memberships)
{
    if (is.list(memberships))
        lapply(memberships, as.gset)
    else
        memberships
}

