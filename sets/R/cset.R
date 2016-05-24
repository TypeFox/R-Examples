#########################
### Customizable Sets ###
#########################

### generator

cset <-
function(gset,
         orderfun = sets_options("orderfun"),
         matchfun = sets_options("matchfun"))
{
    gset <- as.gset(gset)

    ## recreate gset according to user-specified match-fun
    if (!is.null(matchfun)) {
        uni <- !.duplicated_by_matchfun(gset, matchfun)
        gset <-
            .make_gset_by_support_and_memberships(.as.list(.get_support(gset)[uni]),
                                                  .get_memberships(gset)[uni],
                                                  universe = .get_universe(gset),
                                                  bound = .get_bound(gset)
                                                  )
    }

    ## create cset-object
    .make_cset_from_gset_and_orderfun_and_matchfun(gset,
                                                   orderfun,
                                                   matchfun)
}

## convenience function generator for non-vectorized equality predicates
matchfun <-
function(FUN)
    .make_matchfun_from_equalityfun(FUN)

## Disable numeric subscripting (as sets are "unordered" collections of
## elements).  Note that iterating via for() and lapply() still works,
## the former because this [currently, 2007-09-16] directly uses the
## internal list representation and the latter because we provide an
## as.list() method.

`[.cset` <-
function(x, i = x)
{
    ind <- .lookup_elements(x, i, .matchfun(x))
    cset(gset(.as.list(x)[ind], cset_memberships(x)[ind]),
         .orderfun(x),
         .matchfun(x))
}

`[[.cset` <-
function(x, i)
{
    as.set(x)[[i]]
}

`[<-.cset` <-
function(x, i = x, value)
{
    cset(gset(`[<-`(.as.list(x), .lookup_elements(x, i, .matchfun(x)), value),
                memberships = cset_memberships(x)),
         .orderfun(x),
         .matchfun(x))
}

`[[<-.cset` <-
function(x, i, value)
{
    if (!is.character(i) || length(i) > 1L) i <- list(i)
    cset(gset(`[[<-`(.as.list(x), .lookup_elements(x, i, .matchfun(x)), value),
         memberships = cset_memberships(x)),
         .orderfun(x),
         .matchfun(x))
}

### na.omit

na.omit.cset <-
function(object, ...)
{
    m <- .get_memberships(object)
    if (!is.list(m))
        m[is.na(m)] <- 0
    else
        m <- lapply(m, function(i) {
            m2 <- .get_memberships(i)
            m2[is.na(m2) | is.na(i)] <- 0
            .make_gset_from_support_and_memberships(i, m2)
        })
    .make_gset_from_support_and_memberships(object, m)
}

### Ops-method

Ops.cset <-
function(e1, e2)
{
    if(nargs() == 1L) {
        ## dispatch manually for subclasses
        if (inherits(e1, "set"))
            Ops_set(e1, e2, .Generic = .Generic, .Class = .Class)
        if (inherits(e1, "gset"))
            Ops_gset(e1, e2, .Generic = .Generic, .Class = .Class)

        if(!(as.character(.Generic) %in% "!"))
            stop(gettextf("Unary '%s' not defined for \"%s\" objects.",
                          .Generic, .Class),
                 domain = NA)
        return(cset_complement(e1))
    }

    ## dispatch manually for subclasses
    if (inherits(e1, "set") && inherits(e2, "set"))
        Ops_set(e1, e2, .Generic = .Generic, .Class = .Class)
    if (inherits(e1, "gset") && inherits(e2, "gset"))
        Ops_gset(e1, e2, .Generic = .Generic, .Class = .Class)

    if(!(as.character(.Generic)
         %in% c("<", "<=", ">", ">=", "==", "!=",
                "&", "|", "*", "+", "-", "^")))
        stop(gettextf("Generic '%s' not defined for \"%s\" objects.",
                      .Generic, .Class),
             domain = NA)

    if(as.character(.Generic) == "^") {
        if(is.cset(e1) &&
            ((trunc(e2) != e2) || (e2 < 1L)))
            stop("Cartesian product only defined for positive integers.")
        if(is.cset(e2) && (e1 != 2L))
            stop("Operator not defined.")
    }

    switch(.Generic,
           "+"  = cset_sum(e1, e2),
           "|"  = cset_union(e1, e2),
           "-"  = cset_difference(e1, e2),
           "&"  = cset_intersection(e1, e2),
           "*"  = cset_cartesian(e1, e2),
           "<"  = cset_is_proper_subset(e1, e2),
           "<=" = cset_is_subset(e1, e2),
           ">"  = cset_is_proper_subset(e2, e1),
           ">=" = cset_is_subset(e2, e1),
           "==" = cset_is_equal(e1, e2),
           "!=" = !cset_is_equal(e1, e2),
           "^"  = {
               if(is.cset(e2))
                   cset_power(e2)
               else
                   do.call(cset_cartesian, rep.int(list(e1), e2))}
           )

}

### print methods

print.cset <-
function(x, ...)
    print.gset(x, ...)

print.summary.cset <-
function(x, ...)
    print.summary.gset(x, ...)

### format method

format.cset <-
function(x, ...) {
    FUN <- cset_orderfun(x)
    x <- if (isTRUE(gset_is_set(x)))
        .as.list(x)
    else
        .make_list_of_elements_from_cset(x)
    if (is.function(FUN))
        x <- x[FUN(x)]
    else if(is.integer(FUN) && (length(x) == length(FUN)))
        x <- x[FUN]
    .format_set_or_tuple(x, "{", "}", ...)
}

## summary method

summary.cset <-
function(object, ...)
{
    len <- cset_cardinality(object)
    if (na <- is.na(len))
        len <- length.set(object)
    out <- if (len == 0L)
        gettext("The empty set.")
    else if (len == 1L)
        gettext("A customizable set with 1 element.")
    else if (na)
        gettextf("A customizable set with %g elements.", len)
    else
        gettextf("A customizable set with cardinality %g.", len)
    if(!is.null(attr(object, "matchfun")) && !is.null(attr(object, "orderfun")))
        out <- paste(out, "The match and order functions are user-defined.")
    else if(!is.null(attr(object, "matchfun")))
        out <- paste(out, "The match function is user-defined.")
    else if(!is.null(attr(object, "orderfun")))
        out <- paste(out, "The order function is user-defined.")

    .structure(out, class = "summary.cset")
}

### internal stuff

.make_cset_from_gset_and_orderfun_and_matchfun <-
function(gset, orderfun = NULL, matchfun = NULL)
{
    ## make sure that default orderfun and default matchfun are never stored
    if (identical(orderfun, .list_order))
        orderfun <- NULL
    if (identical(matchfun, .exact_match))
        matchfun <- NULL

    ## promote to gset, if only default-funs are specified
    if (is.null(matchfun) && is.null(orderfun))
        return(gset)

    ## create structure (including overwriting gset-class)
    .structure(gset,
               orderfun = orderfun,
               matchfun = matchfun,
               class = "cset")
}

.duplicated_by_matchfun <-
function(x, matchfun)
    duplicated(.as.list(x)[matchfun(x, x)])

.list_unique_by_matchfun <-
function(x, matchfun)
    .as.list(x)[!.duplicated_by_matchfun(x, matchfun)]

.check_matchfun <-
function(l)
{
    matchfun <- cset_matchfun(l[[1]])
    if (!all(sapply(l, function(i) identical(cset_matchfun(i), matchfun))))
        stop("Need same match functions (or none) for all elements.")
    matchfun
}

.check_orderfun <-
function(l)
{
    orderfun <- cset_orderfun(l[[1]])
    if (!all(sapply(l, function(i) identical(cset_orderfun(i), orderfun))))
        NULL
    else
        orderfun
}
