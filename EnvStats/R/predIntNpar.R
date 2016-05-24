predIntNpar <-
function (x, k = m, m = 1, lpl.rank = ifelse(pi.type == "upper", 
    0, 1), n.plus.one.minus.upl.rank = ifelse(pi.type == "lower", 
    0, 1), lb = -Inf, ub = Inf, pi.type = "two-sided") 
{
    if (!is.vector(x, mode = "numeric")) 
        stop("'x' must be a numeric vector.")
    if (!is.numeric(m) || length(m) > 1 || m != trunc(m) || m < 
        1) 
        stop("'m' must be a positive integer")
    if (!is.numeric(k) || length(k) > 1 || k != trunc(k) || k < 
        1 || k > m) 
        stop("'k' must be a positive integer between 1 and 'm'")
    if (any(length.list(lb, ub) != 1) || lb >= ub) 
        stop(paste("'lb' and 'ub' must be scalars,", "and 'lb' must be strictly less than 'ub'"))
    data.name <- deparse(substitute(x))
    if ((bad.obs <- sum(!(x.ok <- is.finite(x)))) > 0) {
        x <- x[x.ok]
        warning(paste(bad.obs, "observations with NA/NaN/Inf in 'x' removed."))
    }
    x <- sort(x)
    n <- length(x)
    if (n < 2 || length(unique(x)) < 2) 
        stop("'x' must contain at least 2 non-missing distinct values.")
    pi.type <- match.arg(pi.type, c("two-sided", "lower", "upper"))
    if (pi.type == "upper") 
        lpl.rank <- 0
    else if (pi.type == "lower") 
        n.plus.one.minus.upl.rank <- 0
    if (!is.numeric(lpl.rank) || length(lpl.rank) > 1 || lpl.rank != 
        trunc(lpl.rank) || lpl.rank < 0 || lpl.rank >= n) 
        stop(paste("'lpl.rank' must be a non-negative integer less than", 
            "the number of non-missing finite values in 'x'"))
    if (pi.type %in% c("two.sided", "lower") & lpl.rank < 1) 
        stop(paste("When pi.type='two.sided' or pi.type='lower',", 
            "'lpl.rank' must be a positive integer"))
    if (!is.numeric(n.plus.one.minus.upl.rank) || length(n.plus.one.minus.upl.rank) > 
        1 || n.plus.one.minus.upl.rank != trunc(n.plus.one.minus.upl.rank) || 
        n.plus.one.minus.upl.rank < 0 || n.plus.one.minus.upl.rank >= 
        n) 
        stop(paste("'n.plus.one.minus.upl.rank' must be a non-negative integer", 
            "less than the number of non-missing finite values in 'x'"))
    if (pi.type %in% c("two.sided", "upper") & n.plus.one.minus.upl.rank < 
        1) 
        stop(paste("When pi.type='two.sided' or pi.type='upper',", 
            "'n.plus.one.minus.upl.rank' must be a positive integer"))
    upl.rank <- n + 1 - n.plus.one.minus.upl.rank
    if (pi.type == "two.sided" && lpl.rank >= upl.rank) 
        stop(paste("Illegal values for 'lpl.rank' and 'n.plus.one.minus.upl.rank'. ", 
            "Make one or both of them smaller."))
    switch(pi.type, `two-sided` = {
        lpl <- x[lpl.rank]
        upl <- x[upl.rank]
        lrl <- lpl.rank
        lru <- upl.rank
    }, lower = {
        lpl <- x[lpl.rank]
        upl <- ub
        lrl <- lpl.rank
        lru <- NULL
    }, upper = {
        lpl <- lb
        upl <- x[upl.rank]
        lrl <- NULL
        lru <- upl.rank
    })
    ret.list <- list(distribution = "None", sample.size = n, 
        data.name = data.name, bad.obs = bad.obs)
    limits <- c(lpl, upl)
    names(limits) <- c("LPL", "UPL")
    vec <- k:m
    conf.level <- sum(exp(lchoose(m - vec + lpl.rank + n.plus.one.minus.upl.rank - 
        1, m - vec) + lchoose(vec + n - lpl.rank - n.plus.one.minus.upl.rank, 
        vec) - lchoose(n + m, m)))
    pi.obj <- list(name = "Prediction", limit.ranks = c(lrl, 
        lru), limits = limits, type = pi.type, method = "Exact", 
        conf.level = conf.level, sample.size = n, k = k, m = m, 
        n.mean = 1)
    oldClass(pi.obj) <- "intervalEstimate"
    ret.list <- c(ret.list, list(interval = pi.obj))
    oldClass(ret.list) <- "estimate"
    ret.list
}
