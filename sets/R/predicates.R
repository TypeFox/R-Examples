## set predicates

is.set <-
function(x)
    inherits(x, "set")

set_is_empty <-
function(x)
{
    if(is.cset(x))
        set_cardinality(x) < 1L
    else
        sapply(x, set_cardinality) < 1L
}

set_is_subset <-
function(x, y)
{
    .help <- function(a, b) set_is_empty(set_complement(b, a))
    if(is.cset(x))
        .help(x, y)
    else
        Vectorize(.help)(x, y)
}

set_is_proper_subset <-
function(x, y)
{
    set_is_subset(x, y) &
    if(is.cset(x))
        length(x) != length(y)
    else
        sapply(x, length) != sapply(y, length)
}

set_is_equal <-
function(x, y)
{
    .help <- function(a, b)
        identical(a, b) ||
        ((length(a) == length(b))
         && (length(set_intersection(a,b)) == length(a)))
    if(is.cset(x))
        .help(x, y)
    else
        Vectorize(.help)(x, y)
}

set_contains_element <-
function(x, e)
{
    if(set_is_empty(x))
        return(FALSE)
    if(is.tuple(e) || is.gset(e) || is_element(e))
        e <- list(e)
    !is.na(.exact_match(e, x))
}

## gset-predicates

is_element <-
function(e)
    inherits(e, "element")

is.gset <-
function(x)
    inherits(x, c("gset", "set"))

gset_is_empty <-
function(x, na.rm = FALSE)
{
    if(is.cset(x))
        gset_cardinality(x, na.rm = na.rm) == 0
    else
        sapply(x, gset_cardinality, na.rm = na.rm) == 0
}

gset_is_subset <-
function(x, y, na.rm = FALSE)
{
    .help <-
        function(a, b) set_is_subset(a, b) &&
    all(unlist(.apply_connector_to_list_of_gsets_and_support(list(a, b),
                                                             .get_support(a),
                                                             `<=`)),
        na.rm = na.rm)
    if(is.cset(x))
        .help(x, y)
    else
        Vectorize(.help)(x, y)
}

gset_is_proper_subset <-
function(x, y, na.rm = FALSE)
{
    gset_is_subset(x, y, na.rm = na.rm) &
    if(is.cset(x))
        cset_cardinality(x, na.rm = na.rm) != cset_cardinality(y, na.rm = na.rm)
    else
        sapply(x, cset_cardinality, na.rm = na.rm) != sapply(y, cset_cardinality, na.rm = na.rm)
}

gset_is_equal <-
function(x, y, na.rm = FALSE)
{
    .help <- function(a, b)
        identical(a, b) ||
        cset_cardinality(a, na.rm = na.rm) ==
            cset_cardinality(b, na.rm = na.rm) &&
    cset_cardinality(gset_intersection(a, b), na.rm = na.rm) ==
        cset_cardinality(a, na.rm = na.rm)
    if(is.cset(x))
        .help(x, y)
    else
        Vectorize(.help)(x, y)
}

gset_contains_element <-
function(x, e)
    set_contains_element(as.set(.make_list_of_elements_from_cset(x)), e(e))

gset_is_multiset <-
function(x, na.rm = FALSE)
{
    m <- .get_memberships(x)
    !is.list(m) && (na.rm && any(m > 1, na.rm = TRUE) ||
                    !na.rm && any(m > 1) && !any(m < 1, na.rm = TRUE))
}

gset_is_crisp <-
gset_is_set_or_multiset <-
function(x, na.rm = FALSE)
{
    m <- .get_memberships(x)
    !is.list(m) && (na.rm && all(m >= 1, na.rm = TRUE) ||
                    !na.rm && all(m >= 1) && !any(m < 1, na.rm = TRUE))
}

gset_is_fuzzy_set <-
function(x, na.rm = FALSE)
{
    m <- .get_memberships(x)
    !is.list(m) && (na.rm && any(m < 1, na.rm = TRUE) ||
                    !na.rm && any(m < 1) && !any(m > 1, na.rm = TRUE))
}

gset_is_set_or_fuzzy_set <-
function(x, na.rm = FALSE)
{
    m <- .get_memberships(x)
    !is.list(m) && (na.rm && all(m <= 1, na.rm = TRUE) ||
                    !na.rm && all(m <= 1) && !any(m > 1, na.rm = TRUE))
}

gset_is_fuzzy_multiset <-
function(x)
    is.list(.get_memberships(x))

gset_is_set <-
function(x, na.rm = FALSE)
{
    m <- .get_memberships(x)
    !is.list(m) && all(m == 1, na.rm = na.rm)
}

gset_has_missings <-
function(x)
{
    M <- .get_memberships(x)
    if (is.list(M))
        any(is.na(unlist(M))) || any(is.na(unlist(lapply(M, .get_memberships))))
    else
        any(is.na(M))
}

## cset-predicates

is.cset <-
function(x)
    inherits(x, c("cset", "gset", "set"))

cset_is_empty <-
function(x, na.rm = FALSE)
{
    if(is.cset(x))
        cset_cardinality(x, na.rm = na.rm) == 0
    else
        sapply(x, cset_cardinality, na.rm = na.rm) == 0
}

cset_is_subset <-
function(x, y, na.rm = FALSE)
{
    .help <-
        function(a, b) set_is_subset(a, b) &&
    all(unlist(.apply_connector_to_list_of_gsets_and_support(list(a, b),
                                                             .get_support(a),
                                                             `<=`)),
        na.rm = na.rm
        )
    if(is.cset(x))
        .help(x, y)
    else
        Vectorize(.help)(x, y)
}

cset_is_proper_subset <-
function(x, y, na.rm = FALSE)
{
    cset_is_subset(x, y, na.rm = na.rm) &
    if(is.cset(x))
        cset_cardinality(x, na.rm = na.rm) !=
            cset_cardinality(y, na.rm = na.rm)
    else
        sapply(x, cset_cardinality, na.rm = na.rm) !=
            sapply(y, cset_cardinality, na.rm = na.rm)
}

cset_is_equal <-
function(x, y, na.rm = FALSE)
{
    .help <- function(a, b)
        identical(a, b) ||
        cset_cardinality(a, na.rm = na.rm) ==
            cset_cardinality(b, na.rm = na.rm) &&
    cset_cardinality(cset_intersection(a, b), na.rm = na.rm) ==
        cset_cardinality(a, na.rm = na.rm)
    if(is.cset(x))
        .help(x, y)
    else
        Vectorize(.help)(x, y)
}

cset_contains_element <-
function(x, e)
{
    if(cset_is_empty(x))
        return(FALSE)
    matchfun <- cset_matchfun(x)
    x <- .make_list_of_elements_from_cset(x)
    e <- e(e)
    if(is.tuple(e) || is.cset(e) || is_element(e))
        e <- list(e)
    ind <- matchfun(e, x)
    if(is.na(ind)) return(FALSE)
    .get_memberships(x[[ind]]) == .get_memberships(e[[1]])
}

cset_is_multiset <-
function(x, na.rm = FALSE)
    gset_is_multiset(x, na.rm = na.rm)

cset_is_crisp <-
cset_is_set_or_multiset <-
function(x, na.rm = FALSE)
    gset_is_crisp(x, na.rm = na.rm)

cset_is_fuzzy_set <-
function(x, na.rm = FALSE)
    gset_is_fuzzy_set(x, na.rm = na.rm)

cset_is_set_or_fuzzy_set <-
function(x, na.rm = FALSE)
    gset_is_set_or_fuzzy_set(x, na.rm = na.rm)

cset_is_fuzzy_multiset <-
function(x)
    gset_is_fuzzy_multiset(x)

cset_is_set <-
function(x, na.rm = FALSE)
    gset_is_set(x, na.rm = na.rm)

cset_has_missings <-
function(x)
    gset_has_missings(x)

## contains_element operator dispatch

"%e%" <-
function(e, x)
{
    if(is.set(x))
        set_contains_element(x, e)
    else if(is.gset(x))
        gset_contains_element(x, e)
    else if(is.cset(x))
        cset_contains_element(x, e)
    else if(is.interval(x))
        interval_contains_element(x, e)
    else
        stop("Predicate undefined.")
}

