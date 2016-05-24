############
### Sets ###
############

### Basic stuff (constructors, print/summary methods)
set <-
function(...) {
    ret <- .set(...)
    n <- names(ret)
    if (!is.null(n) && any(duplicated(n[n != ""])))
        stop("Labels must be unique.")
    ret
}

.set <- function(...)
    .make_set_from_list(.list_sort(.list_unique(list(...))))

print.set <-
function(x, ...)
{
    writeLines(strwrap(format(x, ...), exdent = 1L))
    invisible(x)
}

summary.set <-
function(object, ...)
{
    len <- length(object)
    out <- if (len == 0L)
        gettext("The empty set.")
    else if (len == 1L)
        gettext("A set with 1 element.")
    else
        gettextf("A set with %d elements.", len)
    .structure(out, class = "summary.set")
}

print.summary.set <-
function(x, ...)
{
    writeLines(x)
    invisible(x)
}

format.set <-
function(x, ...) {
    .format_set_or_tuple(x, "{", "}", ...)
}

Math.set <-
function(x, ...)
    as.set(get(.Generic)(unlist(x), ...))


### operators

Ops_set <-
function(e1, e2, .Generic, .Class)
{
    if (missing(e2)) {
        if(!(as.character(.Generic) %in% "!"))
            stop(gettextf("Unary '%s' not defined for \"%s\" objects.",
                          .Generic, .Class),
                 domain = NA)
        return(gset_complement(e1))
    }

    if(!(as.character(.Generic)
         %in% c("<", "<=", ">", ">=", "==", "!=",
                "&", "|", "*", "+", "-", "^")))
        stop(gettextf("Generic '%s' not defined for \"%s\" objects.",
                      .Generic, .Class),
             domain = NA)

    if(as.character(.Generic) == "^") {
        if(is.set(e1) &&
            ((trunc(e2) != e2) || (e2 < 1L)))
            stop("Cartesian product only defined for positive integers.")
        if(is.set(e2) && (e1 != 2L))
            stop("Operator not defined.")
    }

    switch(.Generic,
           "+"  = gset_sum(e1, e2),
           "|"  = set_union(e1, e2),
           "-"  = gset_difference(e1, e2),
           "&"  = set_intersection(e1, e2),
           "*"  = set_cartesian(e1, e2),
           "<"  = set_is_proper_subset(e1, e2),
           "<=" = set_is_subset(e1, e2),
           ">"  = set_is_proper_subset(e2, e1),
           ">=" = set_is_subset(e2, e1),
           "==" = set_is_equal(e1, e2),
           "!=" = !set_is_equal(e1, e2),
           "^"  = {
               if(is.set(e2))
                   set_power(e2)
               else
                   do.call(set_cartesian, rep.int(list(e1), e2))
           }
           )

}

### Subscript methods

## Disable numeric subscripting (as sets are "unordered" collections of
## elements).  Note that iterating via for() and lapply() still works,
## the former because this [currently, 2007-09-16] directly uses the
## internal list representation and the latter because we provide an
## as.list() method.

`[.set` <-
function(x, i = x)
{
    .make_set_from_list(.list_sort(.as.list(x)[.lookup_elements(x, i)]))
}

`[[.set` <-
function(x, i)
{
    if (!is.character(i) || length(i) > 1L || !nzchar(i))
        i <- list(i)
    .as.list(x)[[.lookup_elements(x, i)]]
}

`[<-.set` <-
function(x, i = x, value)
{
    .make_set_from_list(.list_sort(.list_unique(`[<-`(.as.list(x),
                        .lookup_elements(x, i), value))))
}

`[[<-.set` <-
function(x, i, value)
{
    if (!is.character(i) || length(i) > 1L) i <- list(i)
    .make_set_from_list(.list_sort(.list_unique(`[[<-`(.as.list(x),
                        .lookup_elements(x, i), value))))
}

### internal stuff

.make_set_from_list <-
function(x)
   .structure(x, class = c("set", "gset", "cset"))

.format_set_or_tuple <-
function(x, left, right, ...)
{
    nms <- names(x)
    names(x) <- NULL
    SEP <- rep.int("", length(x))
    if (!is.null(nms))
        SEP[nms != ""] <- " = "
    paste(left,
          if (length(x) > 0)
              paste(nms, SEP, LABELS(as.list(x), ...),
                    sep = "", collapse = ", "),
          right,
          sep = "")
}

.set_subset<-
function(x, i)
    as.set(as.list(x)[i])

`.set_replace`<-
function(x, i, value)
{
    ret <- as.list(x)
    ret[i] <- value
    as.set(ret)
}

.set_subset2<-
function(x, i)
    as.list(x)[[i]]

`.set_replace2`<-
function(x, i, value)
{
    ret <- as.list(x)
    ret[[i]] <- value
    as.set(ret)
}

.lookup_elements <-
function(x, values, matchfun = .exact_match)
{
    n <- names(x)
    ind <- if (!is.null(n) && is.character(values))
        matchfun(values, n)
    else {
        values <- .as.list(values)
        matchfun(values, x)
    }
    ind[!is.na(ind)]
}

