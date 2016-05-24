## A simple class for Unicode characters, mostly for printing ...

as.u_char <-
function(x) {
    if(inherits(x, "u_char"))
        return(x)
    if(inherits(x, "u_char_seq"))
        x <- unlist(unclass(x))
    else if(inherits(x, "u_char_range"))
        x <- unlist(lapply(unclass(x), function(e)
                           if(length(e) == 1L) e else seq(e[1L], e[2L]))
                    )
    if(is.double(x) && (x == as.integer(x)))
        x <- as.integer(x)
    if(is.character(x))
        x <- hex(sub("^U\\+", "", x))
    if(is.integer(x) && all(is.na(x) | ((x >= 0L) & (x <= U_max))))
        return(.structure(x, class = "u_char"))
    stop("'x' cannot be coerced to u_char.")
}

`[.u_char` <-
function(x, i)
{
    cl <- oldClass(x)
    y <- NextMethod("[")
    oldClass(y) <- cl
    y
}

as.character.u_char <-
function(x, ...)
    format.u_char(x)

as.data.frame.u_char <-
    as.data.frame.vector

as.integer.u_char <-
function(x, ...)
    unclass(x)

as.list.u_char <-
function(x, ...)
    lapply(seq_along(x), function(i) x[i])

c.u_char <-
function(..., recursive = FALSE)
{
    as.u_char(unlist(lapply(list(...), as.u_char)))
}

format.u_char <-
function(x, ...)
{
    if(!length(x)) return(character())
    ifelse(is.na(x), "<NA>", sprintf("U+%04X", x))
}

print.u_char <-
function(x, ...)
{
    print(noquote(format(x)), ...)
    invisible(x)
}

rep.u_char <-
function(x, times, ...)
    as.u_char(NextMethod("rep"))

unique.u_char <-
function(x, incomparables = FALSE, ...)
    as.u_char(NextMethod("unique"))

## A simple class for Unicode character ranges.

## Using a Lo/Hi data frame might be more efficient, but seems tricky to
## get right ...

as.u_char_range <-
function(x)
{
    if(inherits(x, "u_char_range"))
        return(x)
    if(inherits(x, "u_char"))
        return(.structure(split(x, seq_along(x)),
                          class = "u_char_range"))
    ## <FIXME>
    ## Add more checking eventually:
    ##   We should really have length one or two, and in the latter case
    ##   elements should be sorted ...
    ## </FIXME>
    if(is.character(x)) {
        return(.structure(lapply(strsplit(x, "..", fixed = TRUE),
                                 as.u_char),
                          class = "u_char_range"))
    }
    stop("'x' cannot be coerced to u_char_range.")
}

`[.u_char_range` <-
function(x, i)
{
    y <- unclass(x)[i]
    class(y) <- class(x)
    y
}

as.character.u_char_range <-
function(x, ...)
    format.u_char_range(x)

as.data.frame.u_char_range <-
    as.data.frame.vector

## as.list.u_char_range <-
## function(x, ...)
##     lapply(seq_along(x), function(i) x[i])

c.u_char_range <-
function(..., recursive = FALSE)
{
    ## Not quite perfect ...
    as.u_char_range(unlist(lapply(list(...),
                                  function(e)
                                  format(as.u_char_range(e)))))
}

format.u_char_range <-
function(x, ...)
{
    .structure(as.character(sapply(unclass(x), paste, collapse = "..")),
               names = names(x))
}
    
print.u_char_range <-
function(x, ...)
{
    print(noquote(format(x)), ...)
    invisible(x)
}

rep.u_char_range <-
function(x, times, ...)
{
    x <- format(x)
    as.u_char_range(NextMethod("rep"))
}

unique.u_char_range <-
function(x, incomparables = FALSE, ...)
{
    x <- format(x)
    as.u_char_range(NextMethod("unique"))
}

unlist.u_char_range <-
function(x, recursive = TRUE, use.names = TRUE)
    as.u_char(x)

## A simple class for Unicode character sequences.

as.u_char_seq <-
function(x, sep = NA_character_)
{
    if(inherits(x, "u_char_seq"))
        return(x)
    if(inherits(x, "u_char_range")) {
        x <- unclass(x)
        ind <- as.logical(sapply(x, length) == 2L)
        x[ind] <- lapply(x[ind],
                         function(e) as.u_char(seq.int(e[1L], e[2L],
                                                       1L)))
    } else if(is.list(x))
        x <- lapply(x, as.u_char)
    else if(is.character(x)) {
        ## If sep is empty (or has length zero), take "as is".
        ## Otherwise, if sep is not given, use "," if strings start with
        ## "<", and " " otherwise.
        sep <- as.character(sep)
        x <- if(!length(sep) || identical(sep, ""))
            lapply(x, .str_to_u_char)
        else {
            if(is.na(sep)) {
                sep <- rep.int(" +", length(x))
                ind <- grepl("^<.*>$", x)
                x[ind] <- sub("^<(.*)>$", "\\1", x[ind])
                sep[ind] <- ", *"
                lapply(strsplit(x, sep), as.u_char)
            }
            lapply(strsplit(x, sep), as.u_char)
        }
    } else
        stop("'x' cannot be coerced to u_char_seq.")        
    .structure(x, class = "u_char_seq")
}

`[.u_char_seq` <-
function(x, i)
{
    y <- unclass(x)[i]
    class(y) <- class(x)
    y
}

as.character.u_char_seq <-
function(x, ...)
    format.u_char_seq(x)

as.data.frame.u_char_seq <-
    as.data.frame.vector

## as.list.u_char_seq <-
## function(x, ...)
##     lapply(seq_along(x), function(i) x[i])

c.u_char_seq <-
function(..., recursive = FALSE)
{
    as.u_char_seq(unlist(lapply(list(...), as.u_char_seq),
                         recursive = FALSE))
}

format.u_char_seq <-
function(x, ...)
{
    .structure(sprintf("<%s>", sapply(unclass(x), paste, collapse = ",")),
               names = names(x))
}
    
print.u_char_seq <-
function(x, ...)
{
    print(noquote(format(x)), ...)
    invisible(x)
}

rep.u_char_seq <-
function(x, times, ...)
    as.u_char_seq(NextMethod("rep"))

unique.u_char_seq <-
function(x, incomparables = FALSE, ...)
    as.u_char_seq(NextMethod("unique"))

unlist.u_char_seq <-
function(x, recursive = TRUE, use.names = TRUE)
    as.u_char(x)

## A generic for counting the number of Unicode characters in Unicode
## character data objects.

n_of_u_chars <-
function(x)
    UseMethod("n_of_u_chars")

n_of_u_chars.u_char <-
function(x)
    length(x)

n_of_u_chars.u_char_range <-
function(x)
{
    if(!length(x)) return(integer())
    mat <- sapply(unclass(x), range)
    mat[2L, ] - mat[1L, ] + 1L
}

n_of_u_chars.u_char_seq <-
function(x)
    as.integer(sapply(unclass(x), length))

## A helper function for turning hex codes to integers.

hex <- if(getRversion() >= "2.11.0") {
    function(x) {
        y <- rep.int(NA_integer_, length(x))
        ind <- nzchar(x)
        y[ind] <- strtoi(x[ind], 16L)
        y
    }
} else {
    function(x) {
        y <- rep.int(NA_integer_, length(x))
        re <- "^(0x)?([0123456789abcdeABCDEF]+)$"
        ind <- grepl(re, x)
        y[ind] <-
            as.integer(sapply(parse(text = sub(re, "0x\\2", x[ind])),
                              eval))
        
        y
    }
}

## Unicode codespace ranges from U+0000 to U+10FFFF.
U_max <- hex("10FFFF")

## Matching

u_char_match <-
function(x, table, nomatch = NA_integer_)
{
    ## Match u_char x against a u_char_range table.
    x <- as.u_char(x)
    table <- unclass(as.u_char_range(table))
    y <- rep.int(nomatch, length(x))
    ## <FIXME>
    ## Maybe eventually "improve" by handling the length one cases
    ## directly, but even
    lens <- sapply(table, length)
    ## takes a lot of time.
    ## Remember that we decided not to use a lo/hi matrix representation
    ## for u_char_range objects because getting this into a data frame
    ## did not work ...
    hpos <- cumsum(lens)
    lpos <- 1L + c(0L, hpos[-length(hpos)])
    cps <- unlist(table)
    lo <- cps[lpos]
    hi <- cps[hpos]
    pos <- which(outer(x, lo, ">=") & outer(x, hi, "<="),
                 arr.ind = TRUE)
    ind <- seq_len(nrow(pos))
    y[pos[ind, 1L]] <- pos[ind, 2L]
    ## </FIXME>
    y
}

`%uin%` <-
function(x, table)
    u_char_match(x, table, nomatch = 0L) > 0L

## Internal stuff.

.structure <-
function(x, ...)
    `attributes<-`(x, c(attributes(x), list(...)))
