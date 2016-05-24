############################################################################################
## package 'secr'
## head.R
## first and last n rows of various secr objects
## 2011-10-28
############################################################################################

head.mask <- function (x, n=6L, ...) {
    if (ms(x)) {
        temp <- lapply(x, head, n, ...)
        class(temp) <- class(x)
        temp
    }
    else {
        stopifnot(length(n) == 1L)
        n <- if (n < 0L)
            max(nrow(x) + n, 0L)
        else min(n, nrow(x))
        OK <- seq_len(n)
        subset(x, OK, ...)
    }
}
tail.mask <- function (x, n=6L, ...) {
    if (ms(x)) {
        temp <- lapply(x, tail, n, ...)
        class(temp) <- class(x)
        temp
    }
    else {
        stopifnot(length(n) == 1L)
        nrx <- nrow(x)
        n <- if (n < 0L)
            max(nrx + n, 0L)
        else min(n, nrx)
        OK <- seq.int(to = nrx, length.out = n)
        subset(x, OK, ...)
    }
}

head.Dsurface <- function (x, n=6L, ...) {
    df <- Dsurface.as.data.frame(x, ...)
    head(df, n, ...)
}
tail.Dsurface <- function (x, n=6L, ...) {
    df <- Dsurface.as.data.frame(x, ...)
    tail(df, n, ...)
}

## rely on subset to deal with multi-session traps and capthist
head.traps <- function (x, n=6L, ...) {
    stopifnot(length(n) == 1L)
    n <- if (n < 0L)
        max(nrow(x) + n, 0L)
    else min(n, nrow(x))
    OK <- seq_len(n)
    subset(x, OK, ...)
}
head.capthist <- function (x, n=6L, ...) {
    stopifnot(length(n) == 1L)
    n <- if (n < 0L)
        max(nrow(x) + n, 0L)
    else min(n, nrow(x))
    OK <- seq_len(n)
    subset(x, OK, ...)
}

## must consider varying nrx for multi-session data, so use lapply
tail.traps <- function (x, n=6L, ...) {
    if (ms(x)) {
        temp <- lapply(x, tail, n, ...)
        class(temp) <- class(x)
        temp
    }
    else {
        stopifnot(length(n) == 1L)
        nrx <- nrow(x)
        n <- if (n < 0L)
            max(nrx + n, 0L)
        else min(n, nrx)
        OK <- seq.int(to = nrx, length.out = n)
        subset(x, OK, ...)
    }
}
tail.capthist <- function (x, n=6L, ...) {
    if (ms(x)) {
        temp <- lapply(x, tail, n, ...)
        class(temp) <- class(x)
        temp
    }
    else {
        stopifnot(length(n) == 1L)
        nrx <- nrow(x)
        n <- if (n < 0L)
            max(nrx + n, 0L)
        else min(n, nrx)
        OK <- seq.int(to = nrx, length.out = n)
        subset(x, OK, ...)
    }
}
