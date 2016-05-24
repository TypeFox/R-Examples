##### interval class
### FIXME: actually, multi-intervals -> change class name?

## FIXME: which functions should convert their arguments via
## as.interval?

## constructors

## FIXME: allow for "boxes" (dimensionality > 1)?

## FIXME: what should interval() / integers() etc. return: the domain,
## or an empty set? Currently, intergers/reals/naturals return the domain,
## interval the empty interval.

"%..%" <-
integers <-
function(l = NULL, r = NULL)
{
    if (is.null(l)) l <- -Inf
    if (is.null(r)) r <- Inf
    if (!is.numeric(l) || !is.numeric(r))
        stop("Need numeric bounds.")
    l <- round(l)
    r <- round(r)
    interval(l = l, r = r, domain = "Z")
}

naturals <-
function(l = NULL, r = NULL)
{
    if (is.null(l)) l <- 1
    if (is.null(r)) r <- Inf
    if (!is.numeric(l) || !is.numeric(r))
        stop("Need numeric bounds.")
    l <- round(l)
    r <- round(r)
    if (l < 1L) stop("Positive integers start with 1.")
    interval(l = l, r = r, domain = "N")
}

naturals0 <-
function(l = NULL, r = NULL)
{
    if (is.null(l)) l <- 0
    if (is.null(r)) r <- Inf
    if (!is.numeric(l) || !is.numeric(r))
        stop("Need numeric bounds.")
    l <- round(l)
    r <- round(r)
    if (l < 0) stop("Non-negative integers start with 0.")
    interval(l = l, r = r, domain = "N0")
}

reals <-
function(l = NULL, r = NULL,
         bounds = c("[]", "[)", "(]", "()",
                    "[[", "]]", "][",
                    "open", "closed", "left-open", "right-open",
                    "left-closed", "right-closed")
         )
{
    if (is.null(l)) l <- -Inf
    if (is.null(r)) r <- Inf
    if (!is.numeric(l) || !is.numeric(r))
        stop("Need numeric bounds.")
    interval(l = l, r = r, bounds = bounds, domain = "R")
}

interval <-
function(l = NULL, r = l,
         bounds = c("[]", "[)", "(]", "()",
                    "[[", "]]", "][",
                    "open", "closed", "left-open", "right-open",
                    "left-closed", "right-closed"),
         domain = NULL)
{
    ## bounds
    bounds <- match.arg(bounds)
    bounds <- switch(bounds,
                  "open"=, "]["= "()", "closed"= "[]",
                   "left-open"=, "right-closed"=, "]]"= "(]",
                     "right-open"=, "left-closed"=, "[["= "[)", bounds)

    ## check for empty interval
    if (is.null(l) || is.null(r) || r < l ||
        isTRUE(all.equal(l, r)) && bounds != "[]")
        return(.empty_interval())

    ## domain
    if (is.null(domain)) {
        ## guess domain from mode
        domain <- if ((is.double(l) || is.double(r)) &&
                      !(is.integer(l) && is.infinite(r)) &&
                      !(is.integer(r) && is.infinite(l)))
            "R"
        else if (l > 0)
            "N"
        else if (l > -1L)
            "N0"
        else
            "Z"
    } else {
        ## enforce domain
        domain <- match.arg(domain, c("R", "Z", "N", "N0"))
        if (domain == "R") {
            l <- as.double(l)
            r <- as.double(r)
        } else {
            l <- round(l)
            r <- round(r)
            if (domain == "N" && l < 1L)
                stop("Positive integers start with 1.")
            if (domain == "N0" && l < 0)
                stop("Non-negative integers start with 0.")
        }
    }

    if (!is.numeric(l) || !is.numeric(r))
        stop("Need numeric bounds.")

    ## create interval
    .make_interval_from_list_and_domain(list(list(l = l, r = r,
                                                  lc = substr(bounds, 1L, 1L) == "[",
                                                  rc = substr(bounds, 2L, 2L) == "]"
                                                  )
                                             ),
                                        domain = domain
                                        )
}

## subsetting

`[.interval` <-
function(x, i)
{
    i <- i[i <= length(x)]
    .make_interval_from_list_and_domain(NextMethod(), interval_domain(x))
}

## subscripting
## FIXME: should that be the same than subsetting?

`[[.interval` <-
function(x, i)
    .make_interval_from_list_and_domain(list(NextMethod()), interval_domain(x))

## union = combine

interval_union <-
function(...) c.interval(...)

c.interval <-
function(...)
{
    ## coerce to intervals
    l <- lapply(list(...), as.interval)

    ## remove empty intervals
    l <- l[!sapply(l, interval_is_empty)]
    if (length(l) < 1L) return(.empty_interval())

    ## check domains
    dom <- .maximal_domain(l)

    ## split countable intervals, if any uncountable exists
    countables <- sapply(l, interval_is_countable)
    if (any(countables) && !all(countables))
        l[countables] <- lapply(l[countables], integers2reals)

    l <- unlist(l, recursive = FALSE)

    ## canonicalize
    .make_interval_from_list_and_domain(.canonicalize_interval(l, dom == "R"),
                                        domain = dom)
}

## complement

interval_complement <-
function(x, y = NULL)
{
    x <- as.interval(x)
    x_is_c <- interval_is_countable(x)

    if (!is.null(y)) {
        y <- as.interval(y)
        y_is_c <- interval_is_countable(y)

        ## trivial cases
        if (interval_is_empty(y) || interval_is_empty(x))
            return(y)

        ## promote/reduce x, depending on y
        if (x_is_c && !y_is_c)
            x <- integers2reals(x, min(y), max(y))
        if (!x_is_c && y_is_c)
            x <- reals2integers(x)
        attr(x, "domain") <- attr(y, "domain")

        return(interval_complement(interval_union(x, interval_complement(y))))
    }

    if (interval_is_empty(x))
        return(reals())

    dom <- interval_domain(x)

    i_is_c <- interval_is_countable(x)

    ## bounded parts
    len <- length.interval(x)
    l <- if (len > 1L)
        lapply(seq_len(len - 1L),
               function(i) list(l = max(x[[i]]) + x_is_c,
                                r = min(x[[i + 1L]]) - x_is_c,
                                lc = x_is_c || interval_is_right_open(x[[i]]),
                                rc = x_is_c || interval_is_left_open(x[[i + 1L]])
                                )
               )
    else
        list()

    ## unbounded parts
    if (dom %in% c("R","Z") && interval_is_left_bounded(x) ||
        dom == "N" && min(x) > 1L ||
        dom == "N0" && min(x) > 0L)
        l <- c(list(list(l = if (dom == "N0") 0 else if (dom == "N") 1 else -Inf,
                         r = min(x) - x_is_c,
                         lc = TRUE,
                         rc = x_is_c || interval_is_left_open(x)
                         )
                    ),
               l
               )

    if (interval_is_right_bounded(x))
        l <- c(l,
               list(list(l = max(x) + x_is_c,
                         r = +Inf,
                         rc = TRUE,
                         lc = x_is_c || interval_is_right_open(x)
                         )
                    )
               )

    .make_interval_from_list_and_domain(l, domain = dom)
}

## intersection
interval_intersection <-
function(...)
{
    ## coerce
    len <- length(l <- lapply(list(...), as.interval))

    ## trivial case
    if (any(sapply(l, interval_is_empty)))
        return(.empty_interval())

    ## find minimal domain
    dom <- .minimal_domain(l)
    Dom <- switch(dom, N = naturals(), N0 = naturals0(), Z = integers(), reals())

    ## check for mixture of countables and uncountables
    uc <- sapply(l, interval_is_uncountable)
    if (any(uc) && !all(uc))
        l[uc] <- lapply(l[uc], reals2integers)

    ## workhorse
    IS <- function(x, y)
        interval_complement(interval_union(interval_complement(x),
                                           interval_complement(y),
                                           if (dom == "N") integers(-Inf, 0),
                                           if (dom == "N0") integers(-Inf, -1)
                                           ),
                            Dom)

    ## reduce
    if (len < 1L)
        .empty_interval()
    else if (len < 2L)
        l[[1L]]
    else
        Reduce(IS, l)
}

## symdiff
interval_symdiff <-
function(...)
{
    len <- length(l <- lapply(list(...), as.interval))
    SD <- function(x, y) interval_union(interval_complement(x, y),
                                        interval_complement(y, x))
    if (len < 1L)
        .empty_interval()
    else if (len < 2L)
        l[[1L]]
    else
        Reduce(SD, l)
}

## length
length.interval <-
function(x)
    length(unclass(x))

interval_measure <-
function(x)
    sum(sapply(unclass(x), `[[`, "r") - sapply(unclass(x), `[[`, "l"))

## Summary methods

## FIXME: min/max not defined for open intervals??
## -> we actually compute the infimum  / supremum here.

min.interval <-
function(...)
{
    l <- list(...)
    l <- unclass(do.call(c, l[-length(l)]))
    l[[1L]]$l
}

max.interval <-
function(...)
{
    l <- list(...)
    l <- unclass(do.call(c, l[-length(l)]))
    l[[length(l)]]$r
}

range.interval <- function(...) c(min(...), max(...))

sum.interval <-
function(...)
{
    l <- list(...)
    do.call(interval_sum, l[-length(l)])
}

prod.interval <-
function(...)
{
    l <- list(...)
    do.call(interval_product, l[-length(l)])
}

mean.interval <-
function(...)
{
    l <- list(...)
    do.call(interval_sum, l) / length(l)
}

## operators

Ops.interval <-
function(e1, e2)
{
    if (nargs() == 1L) {
        if(!((ch <- as.character(.Generic)) %in% c("-", "!")))
            stop(gettextf("Unary '%s' not defined for \"%s\" objects.",
                          .Generic, .Class),
                 domain = NA)
        if (ch == "!")
            return(interval_complement(e1))
        else
            return(0L-e1)
    }

    if(!(as.character(.Generic)
         %in% c("<", "<=", ">", ">=", "==", "!=",
                "&", "|", "*", "+", "-", "/", "^")))
        stop(gettextf("Generic '%s' not defined for \"%s\" objects.",
                      .Generic, .Class),
             domain = NA)

    switch(.Generic,
           "|"  = interval_union(e1, e2),
           "&"  = interval_intersection(e1, e2),

           "+"  = interval_sum(e1, e2),
           "-"  = interval_difference(e1, e2),
           "*"  = interval_product(e1, e2),
           "/"  = interval_division(e1, e2),

           "<"  = interval_is_proper_subinterval(e1, e2),
           "<=" = interval_is_subinterval(e1, e2),
           ">"  = interval_is_proper_subinterval(e2, e1),
           ">=" = interval_is_subinterval(e2, e1),
           "==" = interval_is_equal(e1, e2),
           "!=" = !interval_is_equal(e1, e2),

           "^" = if (is.interval(e1) && is.numeric(e2))
               interval_power(e1, e2)
           else if (min(e1) > 0)
               exp(log(e1) * e2)
           else stop("Exponentiation undefined for negative base.")
           )
}

abs.interval <-
function(x)
{
    dom <- interval_domain(x)
    FUN <- function(i) if (i$l >= 0)
        .make_interval_from_list_and_domain(list(i), dom)
    else if (i$r <= 0)
        .make_interval_from_list_and_domain(list(list(l = abs(i$r),
                                                      r = abs(i$l),
                                                      lc = i$rc,
                                                      rc = i$lc
                                                      )
                                                 ),
                                            domain = dom
                                            )
    else
        .make_interval_from_list_and_domain(list(list(l = 0,
                                                      r = max(abs(i$l), i$r),
                                                      lc = TRUE,
                                                      rc = if (abs(i$l) > i$r) i$lc else i$rc
                                                      )
                                                 ),
                                            domain = dom
                                            )

    do.call(c.interval, lapply(unclass(x), FUN))
}

Math.interval <-
function(x, ...)
{
    if (.Generic %in% c("log", "sqrt") && min(x) < 0)
        stop(.Generic, " only defined for positive reals.")

    if (!.Generic %in% c("sqrt", "exp", "log", "log2", "log10",
                         "trunc", "round", "ceiling", "floor",
                         "sin", "cos", "tan", "asin", "acos", "atan"))
        stop(.Generic, " not implemented for intervals.")

    dom <- interval_domain(x)
    FUN <- if (dom == "R")
        function(i) {
            if (.Generic %in% c("sin", "cos", "tan")) {
                if (.Generic %in% c("sin", "cos") &&
                    (i$r - i$l > 2 * pi ||
                     isTRUE(all.equal(i$r - i$l, 2 * pi)))) {
                    l <- -1
                    r <- 1
                } else if (.Generic == "tan" &&
                           (i$r - i$l > pi ||
                            isTRUE(all.equal(i$r - i$l, pi)))) {
                    l <- -Inf
                    r <- Inf
                } else {
                    shift <- ceiling(i$l / (2 * pi)) * 2 * pi
                    l <- i$l - shift
                    r <- i$r - shift

                    p <- pi * -3:3 / 2
                    fp <- p[p > l & p < r |
                            sapply(p, function(i) isTRUE(all.equal(i, l))) |
                            sapply(p, function(i) isTRUE(all.equal(i, r)))]

                    x <- get(.Generic)(c(l, fp, r))
                    l <- min(x)
                    r <- max(x)
                }
            } else if (.Generic == "acos") { ## cos is monotonic decreasing
                l <- get(.Generic)(i$r, ...)
                r <- get(.Generic)(i$l, ...)
            } else { ## all monotonic increasing funs
                l <- get(.Generic)(i$l, ...)
                r <- get(.Generic)(i$r, ...)
            }

            ret <- .make_interval_from_list_and_domain(list(list(l = l,
                                                                 r = r,
                                                                 lc = i$lc,
                                                                 rc = i$rc
                                                                 )
                                                            ),
                                                       domain = dom
                                                       )

            ## dicretize for rounding functions
            if (.Generic %in% c("trunc", "round", "floor")) {
                digits <- list(...)
                if (.Generic != "round" || !length(digits))
                    digits <- 0
                ret <- do.call(c.interval,
                               as.list(seq.int(unclass(ret)[[1]]$l,
                                               unclass(ret)[[1]]$r,
                                               by = 10 ^ -digits[[1]])))
            }
            ret
    } else
        function(i) as.interval(lapply(seq.int(i$l, i$r), .Generic, ...))

    do.call(c.interval, lapply(unclass(x), FUN))
}

## interval arithmetic

.interval_ops <-
function(l, OP)
{
    len <- length(l <- lapply(unclass(l),
                              function(i) integers2reals(as.interval(i))))
    FUN <- function(x, y)
        do.call(c.interval, set_outer(unclass(x), unclass(y), OP))

    if (len < 1L)
        .empty_interval()
    else if (len < 2L)
        l[[1L]]
    else
        Reduce(FUN, unclass(l))
}

.bounds <- function(x, y)
    c("()", "[)", "(]", "[]")[1 + x + 2 * y]

interval_sum <-
function(...)
{
   FUN <- function(x, y)
       interval(x$l + y$l, x$r + y$r,
                .bounds(x$lc && y$lc, x$rc && y$rc))

   .interval_ops(list(...), FUN)
}

interval_difference <-
function(...)
{
    FUN <- function(x, y)
        interval(x$l - y$r, x$r - y$l,
                 .bounds(x$lc && y$rc, x$rc && y$lc))

   .interval_ops(list(...), FUN)

}

interval_product <-
function(...)
{
    FUN <- function(x, y) {
        p <- c(x$l * y$l, x$l * y$r, x$r * y$l, x$r * y$r)
        b <- c(x$lc && y$lc, x$lc && y$rc, x$rc && y$lc, x$rc && y$rc)
        mi <- which.min(p)
        ma <- which.max(p)
        interval(p[mi], p[ma], .bounds(b[mi], b[ma]))
    }

    .interval_ops(list(...), FUN)

}

interval_division <-
function(...)
{
    FUN <- function(x, y) {
        if (0 < y$l || 0 > y$r)
            interval_product(interval(x$l, x$r, .bounds(x$lc, x$rc)),
                             interval(1 / y$r, 1 / y$l, .bounds(y$rc, y$lc)))
        else if (x$l < 0 && 0 < x$r && y$l <= 0 && 0 <= y$r ||
                 y$l == 0 && y$r == 0)
            reals()
        else if (x$r <= 0 && y$l < 0 && y$r == 0)
            interval(x$r / y$l, Inf, .bounds(x$rc && y$lc, TRUE))
        else if (x$r <= 0 && y$l < 0 && 0 < y$r)
            c.interval(interval(-Inf, x$r / y$r, .bounds(TRUE, x$rc && y$rc)),
                       interval(x$r / y$l, Inf, .bounds(x$rc && y$lc, TRUE)))
        else if (x$r <= 0 && y$r > 0 && y$l == 0)
            interval(-Inf, x$r / y$r, .bounds(TRUE, x$rc && y$rc))

        else if (x$l >= 0 && y$l < 0 && y$r == 0)
            interval(-Inf, x$l / y$l, .bounds(TRUE, x$lc && y$lc))
        else if (x$l >= 0 && y$l < 0 && 0 < y$r)
            c.interval(interval(-Inf, x$l / y$l, .bounds(TRUE, x$lc && y$lc)),
                       interval(x$l / y$r, Inf, .bounds(x$lc && y$rc, TRUE)))
        else if (x$l >= 0 && y$r > 0 && y$l == 0)
            interval(x$l / y$r, Inf, .bounds(x$lc && y$rc, TRUE))
        else
            stop("This should not happen.")
    }

    .interval_ops(list(...), FUN)
}

interval_power <-
function(x, n)
{
    x <- as.interval(x)

    if (!isTRUE(all.equal(as.integer(n), n)))
        stop("Power of intervals only defined for integer exponents.")

    if (n == 0)
        return(reals(1,1))

    if (n < 0)
        return(1 / interval_power(x, abs(n)))

    FUN <- function(X)
    {
        if (n %% 2 == 0 && (X$l < 0 || X$r < 0))
            interval_power(abs(x), n)
        else
            .make_interval_from_list_and_domain(list(list(l = X$l ^ n,
                                                          r = X$r ^ n,
                                                          lc = X$lc,
                                                          rc = X$rc
                                                          )
                                                     ),
                                                interval_domain(x)
                                                )
    }

    do.call(c.interval, lapply(unclass(x), FUN))

}

## coercing functions

integers2reals <-
function(x, min = -Inf, max = Inf)
{
    if(interval_is_uncountable(x)) return(x)

    ## split
    x <- unlist(lapply(unclass(x),
                       function(i) unlist(lapply(seq.int(from = max(i$l, min),
                                                         to = min(i$r, max)),
                                                 interval),
                                          recursive = FALSE)
                       ),
                recursive = FALSE)

    ## canonicalize
    .make_interval_from_list_and_domain(.canonicalize_interval(x, reals = TRUE),
                                        "R")
}

reals2integers <-
function(x)
{
    if(interval_is_countable(x)) return(x)
    x <- unclass(x)
    for(i in seq_along(x)) {
        X <- x[[i]]
        if (!X$lc || X$l != round(X$l)) X$l <- floor(X$l) + 1
        if (!X$rc || X$r != round(X$r)) X$r <- trunc(X$r)
        if (X$r < X$l)
            X <- NULL
        else
            X$lc <- X$rc <- TRUE
        x[[i]] <- X
    }

    ## canonicalize
    .make_interval_from_list_and_domain(.canonicalize_interval(x, reals = FALSE),
                                        domain = "Z")
}

as.character.interval <-
function(x, ...)
{
    if (length.interval(x) < 1L)
        return("{}")

    if (interval_is_uncountable(x)) {
        bounds <- if (sets_options("openbounds") == "][")
            c("]", "[", "[", "]")
        else
            c("(", "[", ")", "]")

        ## merge adjacent degenerate intervals into sets
        .merge = function(x, y) {
            i = length(x)
            X <- x[[i]]
            Y <- y[[1]]
            if ((is.set(X) || isTRUE(all.equal(X$l, X$r))) &&
                isTRUE(all.equal(Y$l, Y$r))) {
                if (!is.set(X)) X <- set(X$l)
                c(x[-i], list(c(X, Y$l)))
            } else
            c(x, y)
        }
        l <- Reduce(.merge, lapply(unclass(x), list))

        ## format individual interval
        .format <- function(x) {
            if (is.set(x))
                format(x)
            else if (isTRUE(all.equal(x$l, x$r)))
                paste("{", as.character(x$l), "}", sep = "")
            else
                paste(bounds[1 + x$lc], x$l, ", ", x$r,
                      bounds[3 + x$rc], sep = "")
        }

        paste(sapply(unclass(l), .format), collapse = " U ")

    } else {
        .format <- function(x)
            if (x$l == x$r) x$l else paste(x$l, "..", x$r, sep = "")
        paste(sapply(unclass(x), .format), collapse = ", ")
    }
}

as.interval <-
function(x)
    UseMethod("as.interval")

as.interval.default <- function(x) stop("Not implemented.")
as.interval.numeric <-
function(x)
{
    if (all(is.finite(x)) && isTRUE(all.equal(as.integer(x), x)))
        return(do.call(c.interval, lapply(unclass(x),
                                          function(i) integers(i, i))))
    else
        reals(min(x, na.rm = TRUE), max(x, na.rm = TRUE))
}

as.interval.interval <- function(x) x
as.interval.set <- function(x) as.interval(unclass(x))
as.interval.list <- function(x) do.call(c.interval, x)
as.interval.tuple <- function(x) interval(x[[1]], x[[2]])
as.interval.logical <-
as.interval.NULL <- function(x) .empty_interval()

as.list.interval <- function(x, ...) lapply(seq_along(x), function(i) x[[i]])

make_set_with_order.interval <-
function(x) {
    if (interval_domain(x) == "R") {
        FUN <- function(i)
            if (isTRUE(all.equal(i$l, i$r))) i$l else pair(i$l, i$r)
        .make_set_with_order(.make_set_from_list(lapply(unclass(x), FUN)))
    } else {
        .make_set_with_order(.make_set_from_list(unlist(lapply(unclass(x),
                                                               function(i) seq(i$l, i$r)))))
    }
}

as.double.interval <-
function(x, ...) {
    if (length.interval(x) < 1L)
        return(NULL)
    r <- range(x)
    x <- as.list(reals2integers(x))
    if (r[1] == -Inf) x <- x[-1]
    l <- length(x)
    if (l > 0L && r[2] == Inf) x <- x[-l]
    unlist(lapply(unclass(x), function(i) seq(i$l, i$r)))
}

### print and summary stuff

print.interval <-
function(x, ...)
{
    writeLines(as.character(x))
    invisible(x)
}

summary.interval <-
function(object, ...)
{
    len <- length.interval(object)
    uc <- interval_is_uncountable(object)
    out <- if (len < 1L)
        gettext("The empty interval.")
    else if (len < 2L) {
        if (uc)
            gettext("An interval.")
        else
            gettext("A sequence.")
    } else {
        if (uc)
            gettextf("A union of %d disjoint intervals.", len)
        else
            gettextf("A union of %d sequences.", len)
    }
    .structure(out, class = "summary.interval")
}

print.summary.interval <-
function(x, ...)
{
    writeLines(x)
    invisible(x)
}

format.interval <-
function(x, ...)
    as.character(x)

## extractors
interval_domain <-
function(x)
    attr(x, "domain")

### predicates
interval_is_uncountable <-
function(x)
    interval_domain(x) == "R"

interval_is_countable <- Negate(interval_is_uncountable)

interval_is_bounded <-
interval_is_finite <-
function(x)
    interval_is_left_bounded(x) && interval_is_right_bounded(x)

interval_is_left_bounded <-
function(x)
    is.finite(min(x)) || interval_is_left_open(x)

interval_is_right_bounded <-
function(x)
    is.finite(max(x)) || interval_is_right_open(x)

interval_is_left_unbounded <-
function(x)
    interval_is_left_closed(x) && is.infinite(min(x))

interval_is_right_unbounded <-
function(x)
    interval_is_right_closed(x) && is.infinite(max(x))

interval_is_unbounded <-
function(x)
    interval_is_left_unbounded(x) && interval_is_right_unbounded(x)

interval_is_half_bounded <-
function(x)
    xor(interval_is_left_bounded(x), interval_is_right_bounded(x))

interval_is_degenerate <-
function(x)
    interval_is_empty(x) || isTRUE(all.equal(min(x), max(x)))

interval_is_proper <- Negate(interval_is_degenerate)

interval_is_empty <-
function(x)
    length.interval(x) < 1L

interval_is_closed <-
function(x)
    interval_is_left_closed(x) && interval_is_right_closed(x)

interval_is_left_closed <-
function(x)
    unclass(x)[[1]]$lc

interval_is_right_closed <-
function(x) {
    i <- length.interval(x)
    unclass(x)[[i]]$rc
}

interval_is_left_open <- Negate(interval_is_left_closed)

interval_is_right_open <- Negate(interval_is_right_closed)

interval_is_equal <-
function(x, y)
{
    x <- as.interval(x)
    y <- as.interval(y)
    ctbls <- interval_is_countable(x) && interval_is_countable(y)
    uc <- interval_is_uncountable(x) && interval_is_uncountable(y)
    x <- unclass(x)
    y <- unclass(y)
    length(x) == length(y) && (ctbls || uc) &&
        isTRUE(all.equal(sapply(x, `[[`, "l"), sapply(y, `[[`, "l"))) &&
        isTRUE(all.equal(sapply(x, `[[`, "r"), sapply(y, `[[`, "r"))) &&
        (ctbls || all(sapply(x, `[[`, "lc") == sapply(y, `[[`, "lc")) &&
                  all(sapply(x, `[[`, "rc") == sapply(y, `[[`, "rc")))
}

"%<%" <-
interval_is_less_than <-
function(x, y)
{
    mx <- max(x)
    my <- min(y)
    mx < my ||
      mx == my && (interval_is_right_open(x) || interval_is_left_open(y))
}

"%>%" <-
interval_is_greater_than <-
function(x, y)
{
    mx <- min(x)
    my <- max(y)
    mx > my ||
      mx == my && (interval_is_right_open(x) || interval_is_left_open(y))
}

"%<=%" <-
interval_is_less_than_or_equal <-
function(x, y)
    max(x) <= min(y)

"%>=%" <-
interval_is_greater_than_or_equal <-
function(x, y)
    min(x) >= max(y)

interval_is_proper_subinterval <-
function(x, y)
    interval_is_subinterval(x, y) && !interval_is_equal(x, y)

interval_is_subinterval <-
function(x, y)
{
    x <- as.interval(x)
    y <- as.interval(y)

    ## trivial cases
    if (interval_is_empty(x))
        return(TRUE)

    if (interval_is_empty(y))
        return(FALSE)

    if (interval_is_uncountable(x) &&
        interval_is_proper(x) &&
        interval_is_countable(y))
        return(FALSE)

    ## FIXME: is Real(Inf) == Integer(Inf)?

    ## check if bounds are nested
    if (min(x) < min(y) || max(x) > max(y))
        return(FALSE)

    ## check each element in x against y
    uc <- interval_is_uncountable(x) || interval_is_uncountable(y)
    l <- length(y <- unclass(y))
    j <- 1
    .not_contained <-
        function(a, b) {
            a$l < b$l || a$r > b$r ||
            uc && (isTRUE(all.equal(a$l, b$l)) && a$lc && !b$lc ||
                   isTRUE(all.equal(a$r, b$r)) && a$rc && !b$rc)
        }
    for (i in x) {
        while (.not_contained(i, y[[j]])) {
            j <- j + 1
            if (j > l) return(FALSE)
        }
    }
    TRUE
}

interval_contains_element <-
function(x, y)
    interval_is_subinterval(y, x)

is.interval <- function(x) inherits(x, "interval")

## str method
str.interval <- function(object, ...) str(unclass(object), ...)

### internal functions

.empty_interval <-
function()
    .make_interval_from_list_and_domain(list(), NA)

.make_interval_from_list_and_domain <- function(x, domain = "R")
    .structure(x, class = "interval", domain = domain)

.sort_interval <-
function(x)
    x[order(sapply(unclass(x), `[[`, "l"), sapply(unclass(x), `[[`, "r"))]

.merge_uncountables <-
function(x, y)
{
    len <- length(x)
    if (len < 1L) return(y)
    if (length(y) < 1L) return(x)
    X <- x[[len]]
    Y <- y[[1L]]
    if (Y$l < X$r || isTRUE(all.equal(Y$l, X$r)) && (X$rc || Y$lc)) {
        ## X and Y overlap
        if (isTRUE(all.equal(Y$l, X$l)) && (X$lc != Y$lc))
            ## X left-open, but Y left-closed -> left-close X
            X$lc <- TRUE
        if (Y$r < X$r || isTRUE(all.equal(Y$r, X$r)) && (!Y$rc || X$rc)) {
            ## Y included in X -> discard Y
            c(x[-len], list(X))
        } else {
            ## true overlap -> merge X and Y
            Y$l <- X$l
            Y$lc <- X$lc
            c(x[-len], list(Y))
        }
    } else c(x, y) ## X and Y are disjoint -> return both
 }

.merge_countables <-
function(x, y)
{
    len <- length(x)
    if (len < 1L) return(y)
    if (length(y) < 1L) return(x)
    X <- x[[len]]
    Y <- y[[1L]]
    if (Y$l - X$r < 2L) {
        ## X and Y overlap
        if (Y$r <= X$r)
            ## Y included in X -> discard Y
            c(x[-len], list(X))
        else {
            ## true overlap -> merge X and Y
            Y$l <- X$l
            c(x[-len], list(Y))
        }
    } else c(x, y) ## X and Y are disjoint -> return both
}

.simplify_interval <-
function(x, reals = TRUE)
{
    Reduce(if (reals) .merge_uncountables else .merge_countables,
           lapply(unclass(x), list))
}

.canonicalize_interval <-
function(x, reals = TRUE)
    .simplify_interval(.sort_interval(x), reals = reals)

.minimal_domain <-
function(l)
{
    dom <- rep(sapply(l, interval_domain), sapply(l, length.interval))
    if (any(dom == "N"))
        "N"
    else if (any(dom == "N0"))
        "N0"
    else if (any(dom == "Z"))
        "Z"
    else
        "R"
}

.maximal_domain <-
function(l)
{
    dom <- rep(sapply(l, interval_domain), sapply(l, length.interval))
    if (any(dom == "R"))
        "R"
    else if (any(dom == "Z"))
        "Z"
    else if (any(dom == "N0"))
        "N0"
    else
        "N"
}

