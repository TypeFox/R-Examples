##########################
######### Tuples #########
##########################

### constructors

tuple <-
function(...)
{
    l <- list(...)
    nam <- names(l)
    if (!is.null(nam)) {
        nam <- rep(nam, length.out = length(l))
        names(l) <- nam
    }
    .make_tuple_from_list(l)
}

.make_tuple_from_list <-
function(x)
    .structure(x, class = "tuple")

singleton <-
function(...)
{
    if(length(list(...)) != 1L)
        stop("Exactly one argument expected.")
    tuple(...)
}

pair <-
function(...)
{
    if(length(list(...)) != 2L)
        stop("Exactly two arguments expected.")
    tuple(...)
}

triple <-
function(...)
{
    if(length(list(...)) != 3L)
        stop("Exactly three arguments expected.")
    tuple(...)
}

## quadruple <- function(...) {
##   if(length(list(...)) != 4L)
##     stop("Need four components.")
##   tuple(...)
## }

## quintuple <- function(...) {
##   if(length(list(...)) != 5L)
##     stop("Need five components.")
##   tuple(...)
## }

## sextuple <- function(...) {
##   if(length(list(...)) != 6L)
##     stop("Need six components.")
##   tuple(...)
## }

## septuple <- function(...) {
##   if(length(list(...)) != 7L)
##     stop("Need seven components.")
##   tuple(...)
## }

## octuple <- function(...) {
##   if(length(list(...)) != 8L)
##     stop("Need eight components.")
##   tuple(...)
## }

## nonuple <- function(...) {
##   if(length(list(...)) != 9L)
##     stop("Need nine components.")
##   tuple(...)
## }

## decuple <- function(...) {
##   if(length(list(...)) != 10L)
##     stop("Need ten components.")
##   tuple(...)
## }

### predicates

is.tuple <-
function(x)
    inherits(x, "tuple")

tuple_is_ntuple <-
function(x, n)
    is.tuple(x) && (length(x) == n)

tuple_is_singleton <-
function(x)
    tuple_is_ntuple(x, 1L)

tuple_is_pair <-
function(x)
    tuple_is_ntuple(x, 2L)

tuple_is_triple <-
function(x)
    tuple_is_ntuple(x, 3L)

## tuple_is_quadruple <- function(x)
##   tuple_is_ntuple(x, 4L)

## tuple_is_quintuple <- function(x)
##   tuple_is_ntuple(x, 5L)

## tuple_is_sextuple <- function(x)
##   tuple_is_ntuple(x, 6L)

## tuple_is_septuple  <- function(x)
##   tuple_is_ntuple(x, 7L)

## tuple_is_octuple <- function(x)
##   tuple_is_ntuple(x, 8L)

## tuple_is_nonuple <- function(x)
##   tuple_is_ntuple(x, 9L)

## tuple_is_decuple <- function(x)
##   tuple_is_ntuple(x, 10L)

### operators, methods

c.tuple <-
function(...)
    do.call(tuple,
            unlist(lapply(list(...),
                          function(i) if (is.tuple(i)) i else list(i)),
                   recursive = FALSE)
            )

"[.tuple" <-
function(x, i)
    .make_tuple_from_list(NextMethod("["))

"[<-.tuple" <-
function(x, i, value)
    .make_tuple_from_list(NextMethod("[<-"))

rep.tuple <-
function(x, ...)
    .make_tuple_from_list(rep(.as.list(x), ...))

Ops.tuple <-
function(e1, e2)
{
    if (!missing(e2)) {
        len <- max(length(e1), length(e2))
        e1 <- rep(e1, length.out = len)
        e2 <- rep(e2, length.out = len)
        ret <- lapply(seq_len(len),
                      function(i) do.call(.Generic, list(e1[[i]], e2[[i]])))
    } else {
        ret <- lapply(seq_along(e1),
                      function(i) do.call(.Generic, list(e1[[i]])))
    }

    .make_tuple_from_list(ret)

}

Math.tuple <-
function(x, ...)
    as.tuple(get(.Generic)(unlist(x), ...))

Summary.tuple <-
function(..., na.rm = FALSE)
    do.call(.Generic, c(as.list(...), na.rm = na.rm))

mean.tuple <-
function(x, ...)
{
    x <- as.numeric(x)
    NextMethod()
}

median.tuple <-
function(x, na.rm = FALSE)
{
    x <- as.numeric(x)
    NextMethod()
}

tuple_outer <- set_outer

### output

print.tuple <-
function(x, ...)
{
    writeLines(strwrap(format(x, ...), exdent = 1L))
    invisible(x)
}

format.tuple <-
function(x, ...)
    .format_set_or_tuple(x, "(", ")", ...)

.ntuple <-
function(i, a = FALSE, A = FALSE, plural = FALSE)
{
    if (plural && (a || A))
      stop("Wrong arguments.")
    if (is.tuple(i))
        i <- length(i)
    if (i == 0L)
        paste(if (A) "The" else "the", "empty tuple")
    else if (i > 10L)
        paste(if (A) "A " else if (a) "a " else "",
              i,
              if (plural) "-tuples" else "-tuple",
              sep = "")
    else
        paste(if (a) .TUPLESa[i] else if (A) .TUPLESA[i] else "",
              if (plural) .P_TUPLES[i] else  .S_TUPLES[i], sep = "")
}

### constants

.S_TUPLES <- c("singleton", "pair", "triple", "quadruple", "quintuple",
               "sextuple", "septuple", "octuple", "nonuple", "decuple")
.P_TUPLES <- c("singletons", "pairs", "triples", "quadruples", "quintuples",
               "sextuples", "septuples", "octuples", "nonuples", "decuples")

.TUPLESa <- c("a ", "a ", "a ", "a ", "a ", "a ", "a ", "an ", "a ", "a ")
.TUPLESA <- c("A ", "A ", "A ", "A ", "A ", "A ", "A ", "An ", "A ", "A ")

