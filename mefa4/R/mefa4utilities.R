## utility functions

nonDuplicated <- function(x, y, change.rownames=FALSE, na.rm=FALSE) {
    if (length(dim(x)) != 2)
        stop("'x' must have 2 dimensions")
    z <- deparse(substitute(y))
    if (z %in% colnames(x))
        y <- x[,z]
    if (nrow(x) != length(y))
        stop("non matching arguments 'x' and 'y'")
    keep <- !duplicated(y)
    out <- x[keep,]
    if (na.rm) {
        y <- y[keep]
        out <- out[!is.na(y),]
        y <- y[!is.na(y)]
    } else {
        y <- y[keep]
    }
    if (change.rownames)
        rownames(out) <- y
    out
}

Melt <-
function(x)
{
    if (inherits(x, "Mefa"))
        x <- x@xtab
    if (inherits(x, "mefa"))
        x <- as(x$xtab, "dgCMatrix")
    if (inherits(x, "matrix"))
        x <- as(x, "dgCMatrix")
    if (inherits(x, "sparseMatrix")) {
        if (is.null(rownames(x)))
            stop("missing row names")
        if (is.null(colnames(x)))
            stop("missing column names")
        x <- as(x, "dgTMatrix")
        rows <- x@i + 1L
        cols <- x@j + 1L
        y <- x@x
        out <- data.frame(rows = factor(x@Dimnames[[1]][rows], 
            levels=x@Dimnames[[1]]), 
            cols = factor(x@Dimnames[[2]][cols], 
            levels=x@Dimnames[[2]]), 
            value = y)
    } else if (is.list(x) && all(sapply(x, function(z) 
        inherits(z, "sparseMatrix")))) {
        if (!all(sapply(x[-1], function(z) 
            identical(z@Dimnames, x[[1]]@Dimnames))))
            stop("dimnames of list elements must be identical")
        n <- length(x)
        X <- rows <- cols <- y <- vector("list", n)
        for (k in 1:n) {
            X[[k]] <- as(x[[k]], "dgTMatrix")
            rows[[k]] <- X[[k]]@i + 1L
            cols[[k]] <- X[[k]]@j + 1L
            y[[k]] <- X[[k]]@x
        }
        out <- data.frame(rows = factor(x[[1]]@Dimnames[[1]][unlist(rows)], 
            levels=x[[1]]@Dimnames[[1]]), 
            cols = factor(x[[1]]@Dimnames[[2]][unlist(cols)], 
            levels=x[[1]]@Dimnames[[2]]), 
            segm = as.factor(rep(names(x), sapply(y, length))),
            value = unlist(y))
    } else stop("object class not appropriate")
    out
}

## atatch a time stamp to file names etc.
pasteDate <- 
function(..., sep = " ", collapse = NULL, sep.date = sep)
{
    out <- paste(..., sep = sep, collapse = collapse)
    paste(out, Sys.Date(), sep = sep.date)
}
paste0date <- 
function(..., collapse = NULL)
{
    paste0(paste0(..., collapse = collapse), Sys.Date())
}

## clear up species names
nameAlnum <- 
function(x, capitalize=c("asis", "first", "none", "all", "mixed"), collapse=" ") 
{
    capitalize <- match.arg(capitalize)
    .capwords <- function(x) {
        s <- strsplit(x, " ")[[1]]
        paste(toupper(substring(s, 1, 1)), substring(s, 2),
              sep = "", collapse = " ")
    }
    capwords <- function(x) sapply(x, .capwords)
    if (capitalize == "first") {
        x <- tolower(x)
        capitalize <- "mixed"
    }
    f <- switch(capitalize,
        "asis"=function(x) return(x), 
        "mixed"=capwords, 
        "none"=tolower, 
        "all"=toupper)
    sapply(x, function(z) {
        paste0(f(strsplit(gsub("[^[:alnum:] ]", "", z), " +")[[1]]), collapse=collapse)
    }, USE.NAMES = !is.null(names(x)))
}

## compare sets
compare_sets <- function(x, y) {
    x <- as.factor(x)
    y <- as.factor(y)
    xl <- levels(x)
    yl <- levels(y)
    xa <- levels(droplevels(x))
    ya <- levels(droplevels(y))
    lab <- c(xlength=length(xl), ylength=length(yl),
        intersect=length(intersect(xl, yl)),
        union=length(union(xl, yl)),
        xbutnoty=length(setdiff(xl, yl)),
        ybutnotx=length(setdiff(yl, xl)))
    act <- c(xlength=length(xa), ylength=length(ya),
        intersect=length(intersect(xa, ya)),
        union=length(union(xa, ya)),
        xbutnoty=length(setdiff(xa, ya)),
        ybutnotx=length(setdiff(ya, xa)))
    rbind(labels=lab, unique=act)
}

## find max/min index and value for rows of a matrix
find_max <- function(x) {
    if (is.null(dim(x))) 
        stop("x must be matrix like object with dim attribute")
    if (is.null(colnames(x))) 
        colnames(x) <- paste0("X", seq_len(ncol(x)))
    tmp <- pbapply(x, 1, function(z) {
        i <- which.max(z)
        if (length(i))
            c(i, z[i]) else c(NA, NA)
    })
    i <- factor(tmp[1, ], levels = seq_len(ncol(x)))
    levels(i) <- colnames(x)
    out <- data.frame(index = i, value = tmp[2, ])
    rownames(out) <- rownames(x)
    out
}
find_min <- function(x) {
    out <- find_max(-1 * x)
    out$value <- -1 * out$value
    out
}
