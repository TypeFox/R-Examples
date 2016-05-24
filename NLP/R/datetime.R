ISO_8601_datetime_components <-
    c("year", "mon", "day", "hour", "min", "sec", "tzd")

parse_ISO_8601_datetime <-
function(x)
{
    x <- as.character(x)
    n <- length(x)
    y <- matrix("", n, 7L)
    dimnames(y) <- list(x, ISO_8601_datetime_components)

    pos <- seq_along(x)

    bad <- (is.na(x) |
            (x == "") |
            ((nzchar(x) > 10L) && (substring(x, 11L, 11L) != "T")))
    if(any(bad)) {
        pos <- pos[!bad]
        x <- x[pos]
    }

    dates <- substring(x, 1L, 10L)
    pat <- "^([[:digit:]]{4})(-[[:digit:]]{2})?(-[[:digit:]]{2})?$"
    m <- regmatches(dates, regexec(pat, dates))
    ind <- sapply(m, length) > 0L
    if(!all(ind)) {
        bad[pos[!ind]] <- TRUE
        pos <- pos[ind]
        x <- x[ind]
        m <- m[ind]
    }
    y[pos, 1L : 3L] <-
        do.call(rbind, m)[, 2L : 4L]

    ind <- (nchar(x) > 10L)
    if(any(ind)) {
        if(!all(ind)) {
            pos <- pos[ind]
            x <- x[ind]
        }
        times <- substring(x, 12L)
        pat <- paste("^",
                     "([[:digit:]]{2}):([[:digit:]]{2})",
                     "(:[[:digit:]]{2}([.][[:digit:]]+)?)?",
                     "(Z|[+-][[:digit:]]{2}:[[:digit:]]{2})",
                     "$",
                     sep = "")
        m <- regmatches(times, regexec(pat, times))
        ind <- sapply(m, length) > 0L
        if(!all(ind))
            bad[pos[!ind]] <- TRUE
        y[pos[ind], 4L : 7L] <-
            do.call(rbind, m[ind])[, c(2L, 3L, 4L, 6L)]
    }

    y[, c(2L, 3L, 6L)] <- substring(y[, c(2L, 3L, 6L)], 2L)

    ## Warn about the bad entries.
    if(any(bad)) {
        warning("Invalid entries:",
                paste("\n ", rownames(y)[bad], collapse = " "))
        y[bad, ] <- ""
    }

    ## If we want year to sec as numeric and tzd as character, we need
    ## to do
    ##   y <- as.data.frame(y, stringsAsFactors = FALSE)
    ## and convert variables 1 to 6: note that this would turn empty to
    ## missing ...

    x <- rownames(y)
    w <- which(y != "", arr.ind = TRUE)
    y <- as.data.frame(y, stringsAsFactors = FALSE)
    y[, 1L : 5L] <- lapply(y[, 1L : 5L], as.integer)
    y[[6L]] <- as.numeric(y[[6L]])
    y <- Map(function(u, v) as.list(u[v]),
             split(y, seq_len(n)),
             split(w[, 2L], factor(w[, 1L], seq_len(n))))
    names(y) <- x
    class(y) <- "ISO_8601_datetime"
    y
}

`[.ISO_8601_datetime` <-
function(x, i)
{
    y <- unclass(x)[i]
    class(y) <- class(x)
    y
}

`$.ISO_8601_datetime` <-
function(x, name)
{
    name <- pmatch(name, ISO_8601_datetime_components)
    as.data.frame(x)[[name]]
}

as.matrix.ISO_8601_datetime <-
function(x, ...)
{
    y <- matrix("", length(x), 7L,
                dimnames = list(names(x), ISO_8601_datetime_components))
    nms <- lapply(x, names)
    y[cbind(rep.int(seq_along(x), sapply(nms, length)),
            match(unlist(nms), ISO_8601_datetime_components))] <-
                as.character(unlist(x))
    y
}

as.data.frame.ISO_8601_datetime <-
function(x, row.names = NULL, optional = FALSE, ...)
{
    y <- as.matrix(x)
    y[y == ""] <- NA_character_
    y <- as.data.frame(y, stringsAsFactors = FALSE)
    y[, 1L : 5L] <- lapply(y[, 1L : 5L], as.integer)
    y[[6L]] <- as.numeric(y[[6L]])
    y
}

as.Date.ISO_8601_datetime <-
function(x, ...)
{
    y <- as.matrix(x)
    y[y == ""] <- NA_character_
    as.Date(sprintf("%s-%s-%s", y[, 1L], y[, 2L], y[, 3L]),
            "%Y-%m-%d")
}

as.POSIXct.ISO_8601_datetime <-
function(x, tz = "", ...)
    as.POSIXct(as.POSIXlt(x))

as.POSIXlt.ISO_8601_datetime <-
function(x, tz = "", ...)
{
    y <- as.matrix(x)
    y[y == ""] <- NA_character_
    offsets <- sub(":", "", y[, 7L])
    offsets[offsets == "Z"] <- "+0000"
    y[, 7L] <- offsets
    strptime(do.call(paste, split(y, col(y))),
             "%Y %m %d %H %M %OS %z",
             tz = "UTC")
}

print.ISO_8601_datetime <-
function(x, ...)
{
    y <- as.matrix(x)
    y <- as.data.frame(y, stringsAsFactors = FALSE)
    print(y)
    invisible(x)
}
