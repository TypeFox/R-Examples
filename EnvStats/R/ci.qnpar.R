ci.qnpar <-
function (x, p, lcl.rank = NULL, ucl.rank = NULL, lb = -Inf, 
    ub = Inf, ci.type = c("two-sided", "lower", "upper"), ci.method = c("exact", 
        "normal.approx"), approx.alpha = 0.05, digits = 0) 
{
    x <- sort(x)
    n <- length(x)
    is.null.lcl.rank <- is.null(lcl.rank)
    is.null.ucl.rank <- is.null(ucl.rank)
    if (!is.null.lcl.rank || !is.null.ucl.rank) {
        if (!is.null.lcl.rank && (length(lcl.rank) > 1 || lcl.rank != 
            trunc(lcl.rank) || lcl.rank < 1 || lcl.rank > n)) 
            stop(paste("'lcl.rank' must be an integer between 1 and", 
                n))
        if (!is.null.ucl.rank && (length(ucl.rank) > 1 || ucl.rank != 
            trunc(ucl.rank) || ucl.rank < 1 || ucl.rank > n)) 
            stop(paste("'ucl.rank' must be an integer between 1 and", 
                n))
        if (!is.null.lcl.rank && is.null.ucl.rank) {
            lcl <- x[lcl.rank]
            ucl <- ub
            ci.type <- "lower"
        }
        else if (is.null.lcl.rank && !is.null.ucl.rank) {
            lcl <- lb
            ucl <- x[ucl.rank]
            ci.type <- "upper"
        }
        else {
            if (lcl.rank >= ucl.rank) 
                stop("'lcl.rank' must be stricly less than 'ucl.rank'")
            lcl <- x[lcl.rank]
            ucl <- x[ucl.rank]
            ci.type <- "two-sided"
        }
        ci.method <- "exact"
    }
    else {
        ci.type <- match.arg(ci.type)
        ci.method <- match.arg(ci.method)
        if (ci.method == "exact") {
            switch(ci.type, `two-sided` = {
                ao2 <- approx.alpha/2
                lcl.rank <- qbinom(ao2, n, p) + 1
                ucl.rank <- qbinom(1 - ao2, n, p)
                if ((pbinom(ucl.rank, n, p) - pbinom(lcl.rank - 
                  1, n, p)) <= 1 - approx.alpha) ucl.rank <- min(n, 
                  ucl.rank + 1)
                lcl <- x[lcl.rank]
                ucl <- x[ucl.rank]
            }, lower = {
                lcl.rank <- qbinom(approx.alpha, n, p) + 1
                lcl <- x[lcl.rank]
                ucl <- ub
            }, upper = {
                ucl.rank <- qbinom(1 - approx.alpha, n, p)
                lcl <- lb
                ucl <- x[ucl.rank]
            })
        }
        else {
            vec <- round(ci.normal.approx(n * p, sqrt(n * p * 
                (1 - p)), n, n - 1, ci.type, approx.alpha)$limits, 
                0)
            lcl.rank <- max(1, vec[1])
            ucl.rank <- min(n, vec[2])
            switch(ci.type, `two-sided` = {
                if ((pbinom(ucl.rank, n, p) - pbinom(lcl.rank - 
                  1, n, p)) <= 1 - approx.alpha) ucl.rank <- min(n, 
                  ucl.rank + 1)
                lcl <- x[lcl.rank]
                ucl <- x[ucl.rank]
            }, lower = {
                lcl <- x[lcl.rank]
                ucl <- ub
            }, upper = {
                lcl <- lb
                ucl <- x[ucl.rank]
            })
        }
    }
    conf.level <- switch(ci.type, `two-sided` = pbinom(ucl.rank - 
        1, n, p) - pbinom(lcl.rank - 1, n, p), lower = 1 - pbinom(lcl.rank - 
        1, n, p), upper = pbinom(ucl.rank - 1, n, p))
    ci.limits <- c(lcl, ucl)
    names(ci.limits) <- c("LCL", "UCL")
    pct <- round(100 * p, digits)
    ci.parameter <- paste(pct, number.suffix(pct), " %ile", sep = "")
    ret.obj <- list(name = "Confidence", parameter = ci.parameter, 
        limit.ranks = c(lcl.rank, ucl.rank), limits = ci.limits, 
        type = ci.type, method = ci.method, conf.level = conf.level, 
        sample.size = n)
    oldClass(ret.obj) <- "intervalEstimate"
    ret.obj
}
