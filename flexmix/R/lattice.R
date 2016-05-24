#
#  Copyright (C) Deepayan Sarkar
#  Internal code copied from package lattice for use in flexmix
#

hist.constructor <- 
function (x, breaks, include.lowest = TRUE, right = TRUE, ...) 
{
    if (is.numeric(breaks) && length(breaks) > 1) 
        hist(as.numeric(x), breaks = breaks, plot = FALSE, include.lowest = include.lowest, 
            right = right)
    else hist(as.numeric(x), breaks = breaks, right = right, 
        plot = FALSE)
}

checkArgsAndCall <- 
function (FUN, args) 
{
    if (!("..." %in% names(formals(FUN)))) 
        args <- args[intersect(names(args), names(formals(FUN)))]
    do.call(FUN, args)
}

formattedTicksAndLabels <- 
function (x, at = FALSE, used.at = NULL, labels = FALSE, logsc = FALSE, 
    ..., num.limit = NULL, abbreviate = NULL, minlength = 4, 
    format.posixt = NULL, equispaced.log = TRUE) 
{
    rng <- if (length(x) == 2) 
        as.numeric(x)
    else range(as.numeric(x))
    if (is.logical(logsc) && logsc) 
        logsc <- 10
    have.log <- !is.logical(logsc)
    if (have.log) 
        logbase <- if (is.numeric(logsc)) 
            logsc
        else if (logsc == "e") 
            exp(1)
        else stop("Invalid value of 'log'")
    logpaste <- if (have.log) 
        paste(as.character(logsc), "^", sep = "")
    else ""
    check.overlap <- if (is.logical(at) && is.logical(labels)) 
        TRUE
    else FALSE
    if (is.logical(at)) {
        at <- if (have.log && !equispaced.log) 
            checkArgsAndCall(axisTicks, list(usr = log10(logbase^rng), 
                log = TRUE, axp = NULL, ...))
        else checkArgsAndCall(pretty, list(x = x[is.finite(x)], 
            ...))
    }
    else if (have.log && (length(at) > 0)) {
        if (is.logical(labels)) 
            labels <- as.character(at)
        at <- log(at, base = logbase)
    }
    if (is.logical(labels)) {
        if (have.log && !equispaced.log) {
            labels <- as.character(at)
            at <- log(at, logbase)
        }
        else labels <- paste(logpaste, format(at, trim = TRUE), 
            sep = "")
    }
    list(at = at, labels = labels, check.overlap = check.overlap, 
        num.limit = rng)
}

calculateAxisComponents <- 
function (x, ..., packet.number, packet.list, abbreviate = NULL, 
    minlength = 4) 
{
    if (all(is.na(x))) 
        return(list(at = numeric(0), labels = numeric(0), check.overlap = TRUE, 
            num.limit = c(0, 1)))
    ans <- formattedTicksAndLabels(x, ...)
    rng <- range(ans$num.limit)
    ok <- ans$at >= min(rng) & ans$at <= max(rng)
    ans$at <- ans$at[ok]
    ans$labels <- ans$labels[ok]
    if (is.logical(abbreviate) && abbreviate) 
        ans$labels <- abbreviate(ans$labels, minlength)
    ans
}

extend.limits <- 
function (lim, length = 1, axs = "r", prop = if (axs == "i") 0 else lattice.getOption("axis.padding")$numeric) 
{
    if (all(is.na(lim))) 
        NA_real_
    else if (is.character(lim)) {
        c(1, length(lim)) + c(-1, 1) * if (axs == "i") 
            0.5
        else lattice.getOption("axis.padding")$factor
    }
    else if (length(lim) == 2) {
        if (lim[1] > lim[2]) {
            ccall <- match.call()
            ccall$lim <- rev(lim)
            ans <- eval.parent(ccall)
            return(rev(ans))
        }
        if (!missing(length) && !missing(prop)) 
            stop("'length' and 'prop' cannot both be specified")
        if (length <= 0) 
            stop("'length' must be positive")
        if (!missing(length)) {
            prop <- (as.numeric(length) - as.numeric(diff(lim)))/(2 * 
                as.numeric(diff(lim)))
        }
        if (lim[1] == lim[2]) 
            lim + 0.5 * c(-length, length)
        else {
            d <- diff(as.numeric(lim))
            lim + prop * d * c(-1, 1)
        }
    }
    else {
        print(lim)
        stop("improper length of 'lim'")
    }
}
