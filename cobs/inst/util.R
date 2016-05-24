### Utility functions --- used for testing, etc

zapNsmall <- function(x) if(is.double(x)) zapsmall(x) else x

summaryCobs <- function(x, level = 0.90, ...)
{
    ## Purpose: something like print(summary( cobs.result ))
    ## ----------------------------------------------------------------------
    ## Arguments: x: result of cobs(); level : to be compatible to old alpha=0.1
    ## ----------------------------------------------------------------------
    ## Author: Martin Maechler, Date: 15 Feb 2002, 16:55
    str(lapply(x, zapNsmall), max.level = 1, ...)
    px <- predict(x, interval = "both", level = level)
    print(as.data.frame(px[, c("cb.lo", "ci.lo", "fit", "ci.up", "cb.up")]),...)
    cat("knots :\n"); print(x$knots, ...)
    cat("coef  :\n"); print(x$coef, ...)
    if(!is.null(x$sic)) {
        print(cbind(lambda = x$pp.lambda, SIC = x$sic), ...)
    }
}

cpuTime <- function(expr)
    cat("Time elapsed:", format(system.time(expr)[3]),"\n")

## a cheap version of sfsmisc::rrange():
robrng <- function(x, coef=1.5)
    boxplot.stats(x, coef = coef, do.conf = FALSE, do.out = FALSE)$stats[c(1, 5)]

