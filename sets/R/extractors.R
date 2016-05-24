## extractors

gset_memberships <-
function(x, filter = NULL)
{
    ex <- if(is.null(filter))
        .get_memberships
    else
        function(i) .get_memberships(i[filter])
    if (is.tuple(x)) lapply(x, ex) else ex(x)
}

gset_support <-
function(x)
    as.set(.get_support(x))

gset_core <-
function(x, na.rm = FALSE)
{
    m <- as.double(sapply(.get_memberships(x), max) >= 1)
    if (na.rm) m[is.na(m)] <- FALSE
    .make_gset_from_support_and_memberships(.get_support(x), m)
}

gset_height <-
function(x, na.rm = FALSE)
    max(unlist(.get_memberships(x)), na.rm = na.rm)

gset_peak <-
function(x, na.rm = FALSE)
{
    m <- as.double(sapply(.get_memberships(x), max) ==
                   gset_height(x, na.rm = TRUE))
    if (na.rm) m[is.na(m)] <- FALSE
    .make_gset_from_support_and_memberships(.get_support(x), m)
}

gset_charfun <-
function(x)
{
    if (!is.gset(x))
        stop("Argument 'x' must be a generalized set.")
    ret <- function(e) {
        if (missing(e)) return(.get_memberships(x))
        if (is_element(e)) e <- list(e)
        ret <- .memberships_for_support(x, e)
        if (is.list(ret) && (length(ret) < 2L))
            ret[[1L]]
        else
            ret
    }
    .structure(ret, class = "gset_charfun")
}

gset_universe <-
function(x)
    .get_universe(x)

gset_bound <-
function(x)
    .get_bound(x)

print.gset_charfun <-
function(x, ...)
{
    writeLines("The characteristic function of a generalized set.")
    invisible(x)
}

### * cset extractors (mostly copies of the gset extractors)

cset_memberships <- gset_memberships

cset_support <- gset_support

cset_core <- gset_core

cset_height <- gset_height

cset_peak <- gset_peak

cset_charfun <-
function(x)
{
    if (!is.cset(x))
        stop("Argument 'x' must be a customizable set.")
    matchfun <- cset_matchfun(x)
    ret <- function(e) {
        if (missing(e)) return(.get_memberships(x))
        if (is_element(e)) e <- list(e)
        ret <- .memberships_for_support(x, e, matchfun)
        if (is.list(ret) && (length(ret) < 2L))
            ret[[1L]]
        else
            ret
    }
    .structure(ret, class = "cset_charfun")
}

cset_universe <- gset_universe

cset_bound <- gset_bound

print.cset_charfun <-
function(x, ...)
{
    writeLines("The characteristic function of a customizable set.")
    invisible(x)
}

cset_orderfun <-
function(x)
{
    FUN <- .orderfun(x)
    if (is.null(FUN))
        .list_order
    else
        FUN
}

"cset_orderfun<-" <-
function(x, value) {
    attr(x, "orderfun") <- value
    x
}

cset_matchfun <-
function(x)
{
    FUN <- .matchfun(x)
    if (is.null(FUN))
        .exact_match
    else
        FUN
}

"cset_matchfun<-" <-
function(x, value)
    cset(x, matchfun = value, orderfun = .orderfun(x))

## internal stuff

.get_memberships <-
function(x)
{
    m <- attr(x, "memberships")
    if (is.null(m)) rep.int(1L, length.set(x)) else m
}

.set_memberships <-
function(x, value)
{
    attr(x, "memberships") <- value
    x
}

.get_fuzzy_multi_memberships <-
function(x)
{
    if (!gset_is_fuzzy_multiset(x))
        lapply(.get_memberships(x), function(i) gset(1L, i))
    else
        .get_memberships(x)
}

.get_universe <-
function(x)
{
    u <- attr(x, "universe")
    if (is.null(u))
        u <- sets_options("universe")
    if (!is.null(u))
        as.set(if (is.interval(u)) as.numeric(u) else eval(u))
    else
        as.set(.get_support(x))
}

.set_universe <-
function(x, value)
{
    attr(x, "universe") <- value
    x
}

.get_bound <-
function(x)
{
    b <- attr(x, "bound")
    if (is.null(b))
        b <- sets_options("bound")
    if (!is.null(b))
        b
    else if (set_is_empty(x))
        return(0)
    else
        max(.multiplicities(.get_memberships(x)), na.rm = TRUE)
}

.set_bound <-
function(x, value)
{
    attr(x, "bound") <- value
    x
}

.get_support <-
function(x)
{
    ## simplify, if all elements are of same class and of length 1:
    cl <- unlist(lapply(x, class))
    len <- lapply(x, length)

    ret <- if (all(cl[1] == cl) && all(len == 1L))
        unlist(x, recursive = FALSE)
    else
        x

    if (is.null(ret))
        list()
    else if (is.recursive(ret))
        .as.list(x)
    else
        ret
}

.multiplicities <-
function(x)
{
    if (is.list(x))
        unlist(lapply(x, function(i) sum(.get_memberships(i))))
    else if(any(x > 1, na.rm = TRUE))
        x
    else
        rep(1L, length(x))
}

.orderfun <-
function(x)
    attr(x, "orderfun")

.matchfun <-
function(x)
    attr(x, "matchfun")

