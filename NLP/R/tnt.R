## A simple class for storing tokens and tags ("tagged tokens").

## Conceptually, a single tagged token is a token/tag pair and our
## Tagged_Token objects are sequences (to allow positional access) of
## tagged tokens, i.e., sequences of pairs.
## The implementation actually uses a "pair" (named list of length two)
## of "slots" giving the token and tag sequences.
## Subscripting via [ extracts subsets of tagged tokens.
## Subscripting via $ extracts one slot.

Tagged_Token_slot_names <- c("token", "tag")

Tagged_Token <-
function(token, tag)
{
    token <- as.character(token)
    tag <- as.character(tag)
    if(length(token) != length(tag))
        stop("arguments must have the same length")
    .Tagged_Token_from_args(token, tag)
}

.Tagged_Token_from_args <-
function(token, tag)
{
    x <- list(token, tag)
    names(x) <- Tagged_Token_slot_names
    .Tagged_Token_from_list(x)
}

.Tagged_Token_from_list <-
function(x)
{
    class(x) <- "Tagged_Token"
    x
}

as.Tagged_Token <-
function(x)
    UseMethod("as.Tagged_Token")

as.Tagged_Token.Tagged_Token <- identity

## <FIXME>
## Should this get a '...'?  (And hence the generic, too?)
as.Tagged_Token.TextDocument <-
function(x)
    tagged_words(x)
## </FIXME>

is.Tagged_Token <-
function(x)
    inherits(x, "Tagged_Token")

`[.Tagged_Token` <-
function(x, i)
    .Tagged_Token_from_list(lapply(unclass(x), `[`, i))

## <TODO>
## Implement eventually ...
`[<-.Tagged_Token` <-
function(x, i, value)
    .NotYetImplemented()
## </TODO>

`[[.Tagged_Token` <-
function(x, i)
    .Tagged_Token_from_list(lapply(unclass(x), `[[`, i))

## <TODO>
## Implement eventually ...
`[[<-.Tagged_Token` <-
function(x, i, value)
    .NotYetImplemented()
## </TODO>

## $.Tagged_Token is not really necessary.
`$<-.Tagged_Token` <-
function(x, name, value)
{
    n <- length(x)
    x <- unclass(x)
    if(is.na(pos <- pmatch(name, Tagged_Token_slot_names)))
        stop("invalid element name")
    value <- as.integer(value)
    if(length(value) != n)
        stop("replacement must have the same length as object")
    x[[pos]] <- value
    .Tagged_Token_from_list(x)
}

as.data.frame.Tagged_Token <-
function(x, row.names = NULL, optional = FALSE, ...)
{
    data.frame(token = x$token, tag = x$tag, row.names = row.names)
}

as.list.Tagged_Token <-
function(x,  ...) 
    lapply(seq_along(x), function(i) x[i])

c.Tagged_Token <-
function(..., recursive = FALSE)
{
    args <- lapply(list(...), function(e) unclass(as.Tagged_Token(e)))
    y <- lapply(Tagged_Token_slot_names,
                function(e) unlist(lapply(args, `[[`, e)))
    names(y) <- Tagged_Token_slot_names
    .Tagged_Token_from_list(y)
}

duplicated.Tagged_Token <-
function(x, incomparables = FALSE, ...)
{
    do.call(`&`, lapply(unclass(x), duplicated))
}

format.Tagged_Token <-
function(x, ...)
{
    sprintf("%s/%s", x$token, x$tag)
}

length.Tagged_Token <-
function(x)
    length(x$token)

names.Tagged_Token <-
function(x)
    NULL

## print.Tagged_Token <-
## function(x, ...)
## {
##     print(format(x, ...))
##     invisible(x)
## }

unique.Tagged_Token <-
function(x, incomparables = FALSE, ...)
    x[!duplicated(x)]
