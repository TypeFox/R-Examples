## A simple string class.

String <-
function(x)
{
    y <- enc2utf8(as.character(x)[[1L]])
    class(y) <- "String"
    y
}

as.String <-
function(x)
    UseMethod("as.String")

as.String.String <- identity

as.String.default <-
function(x)
    String(paste(x, collapse = "\n"))

is.String <-
function(x)
    inherits(x, "String")

print.String <-
function(x, ...)
    writeLines(x)

## Provide a '[' method performing slicing (as we cannot provide S3
## methods for substr, and clearly substrings of strings should be
## strings.

## Note that we have no distinction between spans and span arrays (same
## issue as having no distinction between strings and string arrays in
## base R).  Hence, we take spans to always operate in an array context
## (for now: we could add a drop = FALSE argument to have subscripting
## turn character vectors of length one back to strings again).

`[.String` <-
function(x, i, j)
{
    mysubstring <- function(x, i, j) {
        ## substring() recycles to max length of args only when this is
        ## positive.
        if(!length(i))
            character()
        else
            substring(x, i, j)
    }

    if(missing(j)) {
        if(is.Span(i))
            return(mysubstring(x, i$start, i$end))
        if(is.list(i) && length(i) && all(sapply(i, is.Span)))
            return(lapply(i,
                          function(e)
                          mysubstring(x, e$start, e$end)))
    }
    ## Regular slicing operators in a scalar context.
    String(substr(x, i, j))
}

## More methods?
##
## A popular mailing list discussion item is to use a Java style '+'
## operator for concatenating strings (not uniformly liked as the
## corresponding operation is not commutative).

## For now, provide Python-style string repetition.

`*.String` <-
function(e1, e2)
{
    String(paste(rep.int(e1, e2), collapse = ""))
}
