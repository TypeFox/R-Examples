## set Summary methods
Summary.set <-
function(..., na.rm = FALSE)
{
    l <- list(...)
    if (.Generic == "sum")
        return(Reduce(function(i, j) sum(i, as.numeric(j), na.rm = na.rm),
                      l, 0))
    else if (.Generic == "prod")
        return(Reduce(function(i, j) prod(i, as.numeric(j), na.rm = na.rm),
                      l, 1))
    do.call(.Generic, c(na.rm = na.rm, do.call(set_union, l)))
}

mean.set <-
function(x, ...)
{
    mean(as.numeric(x), ...)
}

median.set <-
function(x, na.rm = FALSE)
{
    median(as.numeric(x), na.rm = na.rm)
}

## gset Summary methods

Summary.gset <-
function(..., na.rm = FALSE)
{
    l <- list(...)
    if (any(sapply(l, gset_is_fuzzy_multiset)))
        stop("Operation not defined for fuzzy multisets.")
    l <- lapply(l, function(i) as.set(as.numeric(i) * .get_memberships(i)))
    do.call(.Generic, c(l, na.rm = na.rm))
}

mean.gset <-
function(x, ..., na.rm = FALSE)
{
    if (gset_is_fuzzy_multiset(x))
        stop("Operation not defined for fuzzy multisets.")

    v <- as.numeric(x)
    m <- .get_memberships(x)
    if (na.rm && any(nas <- is.na(m))) {
        v <- v[!nas]
        m <- m[!nas]
    }

    weighted.mean(v, m, na.rm = na.rm)
}

median.gset <-
function(x, na.rm = FALSE)
{
    if (gset_is_fuzzy_multiset(x))
        stop("Operation not defined for fuzzy multisets.")
    x <- if (gset_is_fuzzy_set(x, na.rm = TRUE))
        as.numeric(x) * .get_memberships(x)
    else {
        n <- as.numeric(x)
        m <- .get_memberships(x)
        n[is.na(m)] <- NA
        m[is.na(m)] <- 1
        rep.int(n, times = m)
    }
    median(x, na.rm = na.rm)
}

## cset Summary methods
## FIXME: can we call Summary.gset directly?

Summary.cset <-
function(..., na.rm = FALSE)
{
    l <- list(...)
    if (any(sapply(l, gset_is_fuzzy_multiset)))
        stop("Operation not defined for fuzzy multisets.")
    l <- lapply(l, function(i) as.set(as.numeric(i) * .get_memberships(i)))
    do.call(.Generic, c(l, na.rm = na.rm))
}

mean.cset <-
function(x, ..., na.rm = FALSE)
    mean.gset(x, ..., na.rm = na.rm)

median.cset <-
function(x, na.rm = FALSE)
    median.gset(x, na.rm = na.rm)
