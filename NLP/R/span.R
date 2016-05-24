## A simple span class for storing start and end integer offsets.

## Apache OpenNLP optionally allows storing types in spans: we use 
## Annotation objects (with no ids or features) for this.

## Conceptually, a single span is a start/end pair and our Span objects
## are sequences (to allow positional access) of spans, i.e., sequences
## of pairs.
## The implementation actually uses a "pair" (named list of length two)
## of "slots" giving the start and end value sequences.
## Subscripting via [ extracts subsets of spans.
## Subscripting via $ extracts one slot.

Span_slot_names <- c("start", "end")

Span <-
function(start, end)
{
    start <- as.integer(start)
    end <- as.integer(end)
    if(length(start) != length(end))
        stop("arguments must have the same length")
    .Span_from_args(start, end)
}

.Span_from_args <-
function(start, end)
    .Span_from_list(list(start = start, end = end))

.Span_from_list <-
function(x)
{
    class(x) <- "Span"
    x
}

as.Span <-
function(x)
    UseMethod("as.Span")

as.Span.Span <- identity

as.Span.Annotation <-
function(x)
    .Span_from_args(x$start, x$end)

is.Span <-
function(x)
    inherits(x, "Span")

`[.Span` <-
function(x, i)
    .Span_from_list(lapply(unclass(x), `[`, i))

## <TODO>
## Implement eventually ...
`[<-.Span` <-
function(x, i, value)
    .NotYetImplemented()
## </TODO>

`[[.Span` <-
function(x, i)
    .Span_from_list(lapply(unclass(x), `[[`, i))

## <TODO>
## Implement eventually ...
`[[<-.Span` <-
function(x, i, value)
    .NotYetImplemented()
## </TODO>

## $.Span is not really necessary.
`$<-.Span` <-
function(x, name, value)
{
    n <- length(x)
    x <- unclass(x)
    if(is.na(pos <- pmatch(name, Span_slot_names)))
        stop("invalid element name")
    value <- as.integer(value)
    if(length(value) != n)
        stop("replacement must have the same length as object")
    x[[pos]] <- value
    .Span_from_list(x)
}

Ops.Span <-
function(e1, e2)
{
    ## Allow to add offsets.
    switch(as.character(.Generic),
           "+" = {
               if(is.Span(e1) && is.numeric(e2) && (length(e2) == 1L))
                   return(Span(e1$start + e2, e1$end + e2))
               if(is.Span(e2) && is.numeric(e1) && (length(e1) == 1L))
                   return(Span(e2$start + e1, e2$end + e1))
               stop("Invalid operands.")
           },
           stop(gettextf("'%s' not defined for \"Span\" objects", 
                         .Generic),
                domain = NA))
}

as.data.frame.Span <-
function(x, row.names = NULL, optional = FALSE, ...)
{
    data.frame(start = x$start, end = x$end, row.names = row.names)
}

as.list.Span <-
function(x,  ...) 
    lapply(seq_along(x), function(i) x[i])

c.Span <-
function(..., recursive = FALSE)
{
    args <- lapply(list(...), function(e) unclass(as.Span(e)))
    y <- lapply(Span_slot_names,
                function(e) unlist(lapply(args, `[[`, e)))
    names(y) <- Span_slot_names
    .Span_from_list(y)
}

duplicated.Span <-
function(x, incomparables = FALSE, ...)
{
    do.call(`&`, lapply(unclass(x), duplicated))
}

format.Span <-
function(x, ...)
{
    format(as.data.frame(x), ...)
}

length.Span <-
function(x)
    length(x$start)

names.Span <-
function(x)
    NULL

print.Span <-
function(x, ...)
{
    print.data.frame(format(x), ..., row.names = FALSE)
    invisible(x)
}

unique.Span <-
function(x, incomparables = FALSE, ...)
    x[!duplicated(x)]
