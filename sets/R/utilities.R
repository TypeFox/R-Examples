### * .cartesian_product

.cartesian_product <-
function(x)
{
    ## Cf. expand.grid().
    out <- vector("list", length(x))
    rep_fac <- 1L
    d <- sapply(x, length)
    orep <- prod(d)
    for(i in seq_along(x)) {
        nx <- d[i]
        orep <- orep / nx
        out[[i]] <-
            as.list(x[[i]])[rep.int(rep.int(seq_len(nx),
                                            rep.int(rep_fac, nx)), orep)]
        rep_fac <- rep_fac * nx
    }
    out
}

### * .make_set_of_tuples_from_list_of_lists

.make_set_of_tuples_from_list_of_lists <-
function(x)
{
    ret <- unlist(x, recursive = FALSE)
    dim(ret) <- c(length(x[[1L]]), length(x))
    as.set(apply(ret, 1L, as.tuple))
}

### make sure that list elements are not destroyed during set unions

.make_list_elements <-
function(i)
    if (!is.cset(i)) list(i) else i

### as.character replacing empty strings by "_", as needed by environments
.as.character <-
function(x)
{
    x <- as.character(x)
    ind <- !is.na(x)
    if (any(ind)) {
        y <- x[ind]
        y[!nzchar(y)] <- "_"
        y <- substr(y, start = 1, stop = 256)
        while(any((n <- nchar(y, "byte")) > 256))
            y <- substr(y, 1, 256 - ceiling(n / 2))
        x[ind] <- y
    }
    x
}

## hash creator and accessors
.create_hash <-
function(s)
{
    if(length.set(s) < 1L) return(NULL)
    s <- .as.character(s)
    h <- new.env()
    for (i in seq.int(from = length(s), to = 1))
        h[[ s[[i]] ]] <- i
    h
}

### exact_match

.make_matchfun_from_equalityfun <-
function(equalityfun)
{
    equalityfun <- match.fun(equalityfun)
    function(x, table)
    {
        table <- .as.list(table)
        hash <- if (sets_options("hash")) .create_hash(table) else NULL
        l <- length(table)
        FUN <- function(i) {
            from <- if (is.null(hash))
                1
            else
                hash[[ .as.character(list(i)) ]]
            if (is.null(from) || l < 1L) return(NA)
            for (ind in seq.int(from = from, to = l))
                if (equalityfun(table[[ind]], i)) break
            if (ind < l || equalityfun(table[[ind]], i))
                ind else NA
        }
        ret <- lapply(.as.list(x), FUN)
        if (length(ret) < 1L) integer() else unlist(ret)
    }
}

.exact_match <-
    .make_matchfun_from_equalityfun(identical)

### .list_order/sort/unique

.list_order <-
function(x, decreasing = FALSE, ...) {
    .as.character <- function(x, ...) {
## FIXME: also serialize S4 objects
        rec <- unlist(lapply(x, is.recursive))
        x[rec] <- lapply(x[rec], serialize, NULL)
        unlist(lapply(x, paste, collapse = ""))
    }

    len <- length(x)
    if (len < 1L)
        return(integer(0L))

    l <- unlist(lapply(x, length))
    cl <- paste(lapply(x, class))

    num <- (cl == "integer" | cl == "numeric") & (l == 1L)
    numind <- numeric(len)
    numind[num] <- unlist(x[num])

    ch <- character(len)
    if (!all(num)) {
        ch[!num] <- .as.character(x[!num])
        loc <- ""
        suppressWarnings(if (capabilities("iconv")) {
            loc <<- Sys.getlocale("LC_COLLATE")
            on.exit(Sys.setlocale("LC_COLLATE", loc))
            Sys.setlocale("LC_COLLATE", "C")
            ch <- iconv(ch, to = "UTF-8")
        })
    }

    order(l, ## length
          cl, ## class
          numind, ## for numeric mode: numeric value
          ch, ## default: character representation
          decreasing = decreasing, ...)
}

.list_sort <-
function(x, decreasing = FALSE, ...)
    .as.list(x)[.list_order(x, decreasing = decreasing, ...)]

.list_unique <-
function(x)
    .as.list(x)[!duplicated(x)]

.domain_is_numeric <-
function(x)
    all(unlist(lapply(x, is.numeric)) &
        (unlist(lapply(x, length)) == 1L))

### fast structure replacement

.structure <- function(x, ...)
    `attributes<-`(x, c(attributes(x), list(...)))
