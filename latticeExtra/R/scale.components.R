

## adapted from the Lattice book by Deepayan Sarkar

xscale.components.logpower <- function(lim, ...) {
    ans <- xscale.components.default(lim, ...)
    ans$bottom$labels$labels <- parse(text = ans$bottom$labels$labels)
    ans
}

yscale.components.logpower <- function(lim, ...) {
    ans <- yscale.components.default(lim, ...)
    ans$left$labels$labels <- parse(text = ans$left$labels$labels)
    ans
}

xscale.components.fractions <- function(lim, logsc = FALSE, ...) {
    ans <- xscale.components.default(lim, logsc = logsc, ...)
    ## get 'at' in data coordinates
    if (identical(logsc, TRUE)) logsc <- 10
    if (identical(logsc, "e")) logsc <- exp(1)
    at <- ans$bottom$labels$at
    if (!identical(logsc, FALSE))
        at <- logsc ^ at
    ans$bottom$labels$labels <- MASS::fractions(at)
    ans
}

yscale.components.fractions <- function(lim, logsc = FALSE, ...) {
    ans <- yscale.components.default(lim, logsc = logsc, ...)
    ## get 'at' in data coordinates
    if (identical(logsc, TRUE)) logsc <- 10
    if (identical(logsc, "e")) logsc <- exp(1)
    at <- ans$left$labels$at
    if (!identical(logsc, FALSE))
        at <- logsc ^ at
    ans$left$labels$labels <- MASS::fractions(at)
    ans
}

## compute nice log-ticks.  This is a version from the Lattice book
## that is not very sophisticated.

logTicksOld <- function (lim, loc = c(1, 5)) {
    ii <- floor(log10(range(lim))) + c(-1, 2)
    main <- 10^(ii[1]:ii[2])
    r <- as.numeric(outer(loc, main, "*"))
    r[lim[1] <= r & r <= lim[2]]
}

## A more sophisticated version that uses the same algorithm used in
## traditional graphics, via axisTicks() - new in R 2.14.0

logTicks <- function (lim, loc = NULL) {
    if (is.null(loc)) axisTicks(log10(lim), log=TRUE)
    else logTicksOld(lim, loc)
}

xscale.components.log <- function(lim, logsc = FALSE, at = NULL, loc = NULL, ...) {
    ans <- xscale.components.default(lim = lim, logsc = logsc, at = at, ...)
    if (is.null(at)) return(ans)
    if (identical(logsc, FALSE)) return(ans)
    logbase <- logsc
    if (identical(logbase, TRUE)) logbase <- 10
    if (identical(logbase, "e")) logbase <- exp(1)
    tick.at <- logTicks(logbase^lim, loc = loc)
    ans$bottom$ticks$at <- log(tick.at, logbase)
    ans$bottom$labels$at <- log(tick.at, logbase)
    ans$bottom$labels$labels <- as.character(tick.at)
    ans
}

yscale.components.log <- function(lim, logsc = FALSE, at = NULL, loc = NULL, ...) {
    ans <- yscale.components.default(lim = lim, logsc = logsc, at = at, ...)
    if (is.null(at)) return(ans)
    if (identical(logsc, FALSE)) return(ans)
    logbase <- logsc
    if (identical(logbase, TRUE)) logbase <- 10
    if (identical(logbase, "e")) logbase <- exp(1)
    tick.at <- logTicks(logbase^lim, loc = loc)
    ans$left$ticks$at <- log(tick.at, logbase)
    ans$left$labels$at <- log(tick.at, logbase)
    ans$left$labels$labels <- as.character(tick.at)
    ans
}

xscale.components.log10.3 <- function(lim, logsc = FALSE, at = NULL, ...) {
    xscale.components.log(lim, logsc = logsc, at = at, loc = c(1, 3)) 
}

yscale.components.log10.3 <- function(lim, logsc = FALSE, at = NULL, ...) {
    yscale.components.log(lim, logsc = logsc, at = at, loc = c(1, 3))
}


# major + minor ticks for powers of 10

xscale.components.log10ticks <- function(lim, logsc = FALSE, at = NULL, ...) {
    ans <- xscale.components.default(lim = lim, logsc = logsc, at = at, ...)
    if (is.null(at)) return(ans)
    if (identical(logsc, FALSE)) return(ans)
    logbase <- logsc
    if (identical(logbase, TRUE)) logbase <- 10
    if (identical(logbase, "e")) logbase <- exp(1)
    tick.at <- logTicks(logbase^lim, loc = 1:9)
    tick.at.major <- logTicks(logbase^lim, loc = 1)
    major <- tick.at %in% tick.at.major
    ans$bottom$ticks$at <- log(tick.at, logbase)
    ans$bottom$ticks$tck <- ifelse(major, 1, 0.5)
    ans$bottom$labels$at <- log(tick.at, logbase)
    ans$bottom$labels$labels <- as.character(tick.at)
    ans$bottom$labels$labels[!major] <- ""
    ans$bottom$labels$check.overlap <- FALSE
    ans
}

yscale.components.log10ticks <- function(lim, logsc = FALSE, at = NULL, ...) {
    ans <- yscale.components.default(lim = lim, logsc = logsc, at = at, ...)
    if (is.null(at)) return(ans)
    if (identical(logsc, FALSE)) return(ans)
    logbase <- logsc
    if (identical(logbase, TRUE)) logbase <- 10
    if (identical(logbase, "e")) logbase <- exp(1)
    tick.at <- logTicks(logbase^lim, loc = 1:9)
    tick.at.major <- logTicks(logbase^lim, loc = 1)
    major <- tick.at %in% tick.at.major
    ans$left$ticks$at <- log(tick.at, logbase)
    ans$left$ticks$tck <- ifelse(major, 1, 0.5)
    ans$left$labels$at <- log(tick.at, logbase)
    ans$left$labels$labels <- as.character(tick.at)
    ans$left$labels$labels[!major] <- ""
    ans$left$labels$check.overlap <- FALSE
    ans
}


## major + minor ticks (e.g. for date/time axes):

xscale.components.subticks <-
    function(lim, ..., n = 5, n2 = n * 5, min.n2 = n + 5)
{
    ans <- xscale.components.default(lim = lim, ..., n = n)
    ans2 <- xscale.components.default(lim = lim, ..., n = n2, min.n = min.n2)
    ticks <- ans$bottom$ticks$at
    ticks2 <- ans2$bottom$ticks$at
    ticks2 <- ticks2[!(ticks2 %in% ticks)]
    ans$bottom$ticks$at <- c(ticks, ticks2)
    ans$bottom$ticks$tck <- c(rep(1, length(ticks)),
                              rep(0.5, length(ticks2)))
    ans$bottom$labels$at <- ans$bottom$ticks$at
    ans$bottom$labels$labels <- c(ans$bottom$labels$labels,
                                  rep(" ", length(ticks2)))
    ans$bottom$labels$check.overlap <- FALSE
    ans
}

yscale.components.subticks <-
    function(lim, ..., n = 5, n2 = n * 5, min.n2 = n + 5)
{
    ans <- yscale.components.default(lim = lim, ..., n = n)
    ans2 <- yscale.components.default(lim = lim, ..., n = n2, min.n = min.n2)
    ticks <- ans$left$ticks$at
    ticks2 <- ans2$left$ticks$at
    ticks2 <- ticks2[!(ticks2 %in% ticks)]
    ans$left$ticks$at <- c(ticks, ticks2)
    ans$left$ticks$tck <- c(rep(1, length(ticks)),
                            rep(0.5, length(ticks2)))
    ans$left$labels$at <- ans$left$ticks$at
    ans$left$labels$labels <- c(ans$left$labels$labels,
                                rep(" ", length(ticks2)))
    ans$left$labels$check.overlap <- FALSE
    ans
}
