
`shuffleStrata` <- function(strata, type, mirror = FALSE, start = NULL,
                            flip = NULL, nrow, ncol, start.row = NULL,
                            start.col = NULL) {
    ## drop unused levels
    strata <- droplevels(strata)
    LEVS <- levels(strata)
    lev <- nlevels(strata)
    ngr <- length(strata) / lev
    SEQ <- seq_len(lev)
    sp <- split(out <- seq_along(strata), strata)
    perm <- if(type == "free") {
        shuffleFree(lev, lev)
    } else if (type == "series") {
        shuffleSeries(SEQ, mirror = mirror, start = start,
                      flip = flip)
    } else if (type == "grid") {
        shuffleGrid(nrow = nrow, ncol = ncol, mirror = mirror,
                    start.row = start.row, start.col = start.col,
                    flip = flip)
    } else {
        stop("Invalid permutation type.")
    }
    for(i in SEQ) {
        want <- which(strata == LEVS[i])
        out[want] <- sp[[perm[i]]]
    }
    out
}

`shuffleGrid` <- function(nrow, ncol, mirror = FALSE, start.row = NULL,
                          start.col = NULL, flip = NULL) {
    if(is.null(start.row))
        start.row <- shuffleFree(nrow, 1L)
    if(is.null(start.col))
        start.col <- shuffleFree(ncol, 1L)
    ir <- seq(start.row, length=nrow) %% nrow
    ic <- seq(start.col, length=ncol) %% ncol
    if(!is.null(flip) && mirror) {
        if(any(flip)) {
            if(flip[1L])
                ir <- rev(ir)
            if(flip[2L])
                ic <- rev(ic)
        }
    } else {
        if (mirror) {
            if (runif(1L) < 0.5)
                ir <- rev(ir)
            if (runif(1L) < 0.5)
                ic <- rev(ic)
        }
    }
    rep(ic, each=nrow) * nrow + rep(ir, len=nrow*ncol) + 1L
}

`shuffleSeries` <- function(x, mirror = FALSE, start = NULL,
                            flip = NULL) {
    n <- length(x)
    if(is.null(start))
        start <- shuffleFree(n, 1L)
    out <- seq(start, length = n) %% n + 1L
    if(!is.null(flip) && mirror) {
        if(flip)
            out <- rev(out)
    } else {
        if(mirror && runif(1L) < 0.5)
            out <- rev(out)
    }
    x[out]
}

`shuffleFree` <- function(x, size) {
    sample.int(x, size, replace = FALSE)
}

## wrapper function when shuffling without any strata at all at any level
`shuffleNoStrata` <- function(n, control) {
    type <- control$within$type
    switch(type,
           "free" = shuffleFree(n, n),
           "series" = shuffleSeries(seq_len(n), mirror = control$within$mirror),
           "grid" = shuffleGrid(nrow = control$within$nrow,
           ncol = control$within$ncol, mirror = control$within$mirror),
           "none" = seq_len(n)
           )
}
