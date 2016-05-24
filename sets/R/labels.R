### LABELS
###
## We need a function that produces "nice" labels from any object mainly
## for printing, but also for dimnames.  This function should do more
## than as.character(), and less than format ...
##
## Normally, one uses LABELS() and adds extensions by writing methods
## for LABEL().
##
## What we do in LABELS() is the following:
##
## 1. transform the given object to a list
## 2. check names attribute; if any, use these as default
## 3. for all components with empty name, use LABEL() to compute a
##    "simple" representation
## 4. optionally, truncate strings to specified length
## 5. optionally, apply make.unique() to the result
##
## Generally, LABEL() uses format() if the argument is of length 1, and
## creates a type specification otherwise.
## Exception: we also accept "small" sets and pairs since they can well
## be distinguished even if they are nested.  Currently, "small" means a
## length of 5 which is sort of ad-hoc.

LABELS <-
function(x, max_width = NULL, dots = "...", unique = FALSE, limit = NULL, ...)
{
    x <- as.list(x)
    l <- length(x)

    ## recycle max_width and dots as needed
    if (!is.null(max_width))
        max_width <- rep(max_width, length.out = l)
    dots <- rep(dots, length.out = l)

    ## check existing labels
    ret <- names(x)
    if (is.null(ret))
        ret <- rep.int("", l)

    ## create a label for components without given one
    empty <- is.na(ret) | (ret == "")
    if (any(empty))
        ret[empty] <- sapply(x[empty], LABEL, limit, ...)

    ## check maximum width (max_width == NULL => unbounded)
    if (!is.null(max_width)) {
        too_long <- nchar(ret, "width") > max_width
        if (any(too_long)) {
            ret[too_long] <- strtrim(ret[too_long], max_width[too_long])

            ## possibly add dots
            if (!is.null(dots))
                ret[too_long] <- paste(ret[too_long], dots[too_long], sep = "")
          }
    }

    if (unique)
      ret <- make.unique(ret)

    ret
}

LABEL <-
function(x, limit = NULL, ...)
    UseMethod("LABEL")

LABEL.default <-
function(x, limit = NULL, ...)
    paste("<<", class(x)[1L], ">>", sep = "")

LABEL.matrix <-
function(x, limit = NULL, ...)
    sprintf("<<%ix%i matrix>>", nrow(x), ncol(x))

LABEL.numeric <-
LABEL.factor <-
LABEL.integer <-
LABEL.logical <-
function(x, limit = NULL, ...) {
    if (is.null(limit))
        limit <- 2L
    .format_or_class(x, limit, ...)
}

LABEL.character <-
function(x, limit = NULL, quote = sets_options("quote"), ...) {
    if (is.null(limit))
        limit <- 2L
    if (quote)
        x <- ifelse(is.na(x), x, paste("\"", x, "\"", sep = ""))
    .format_or_class(x, limit, ...)
}

LABEL.list <-
function(x, limit = NULL, ...) {
    if (is.null(limit))
        limit <- 1L
    .format_or_class(x, limit, ...)
}

LABEL.set <-
LABEL.gset <-
LABEL.cset <-
LABEL.tuple <-
LABEL.interval <-
function(x, limit = NULL, ...) {
    if (is.null(limit))
        limit <- 6L
    .format_or_class(x, limit, ...)
}

.format_or_class <-
function(x, limit, ...)
{
    l <- length.set(x)
    if (l < limit) {
        if (is.integer(x))
            format(ifelse(is.na(x), x, paste(x, "L", sep = "")), ...)
        else
            format(x, ...)
    } else
        paste("<<", class(x)[1L], "(", l, ")>>", sep = "")
}
