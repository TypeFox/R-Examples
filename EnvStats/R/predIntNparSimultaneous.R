predIntNparSimultaneous <-
function (x, n.median = 1, k = 1, m = 2, r = 1, rule = "k.of.m", 
    lpl.rank = ifelse(pi.type == "upper", 0, 1), n.plus.one.minus.upl.rank = ifelse(pi.type == 
        "lower", 0, 1), lb = -Inf, ub = Inf, pi.type = "upper", 
    integrate.args.list = NULL) 
{
    if (!is.vector(x, mode = "numeric")) 
        stop("'x' must be a numeric vector.")
    if (!is.vector(n.median, mode = "numeric") || length(n.median) != 
        1 || n.median != trunc(n.median) || n.median < 1 || !is.odd(n.median)) 
        stop("'n.median' must be a positive odd integer")
    rule <- match.arg(rule, c("k.of.m", "CA", "Modified.CA"))
    switch(rule, k.of.m = {
        if (!is.vector(k, mode = "numeric") || length(k) != 1 || 
            k != trunc(k) || k < 1 || !is.vector(m, mode = "numeric") || 
            length(m) != 1 || m != trunc(m) || m < 1 || !is.vector(r, 
            mode = "numeric") || length(r) != 1 || r != trunc(r) || 
            r < 1 || k > m) stop(paste("'k', 'm', and 'r' must be positive integers,", 
            "and 'k' must be between 1 and 'm'"))
    }, CA = {
        if (!is.vector(m, mode = "numeric") || length(m) != 1 || 
            m != trunc(m) || m < 1 || !is.vector(r, mode = "numeric") || 
            length(r) != 1 || r != trunc(r) || r < 1) stop("'m', and 'r' must be positive integers")
        k <- "First.or.all.of.next.m.minus.one"
    }, Modified.CA = {
        if (!is.vector(r, mode = "numeric") || length(r) != 1 || 
            r != trunc(r) || r < 1) stop("'r' must be a positive integer")
        k <- "First.or.at.least.two.of.next.three"
        m <- 4
    })
    if (any(length.list(lb, ub) != 1) || lb >= ub) 
        stop(paste("'lb' and 'ub' must be scalars,", "and 'lb' must be strictly less than 'ub'"))
    data.name <- deparse(substitute(x))
    if ((bad.obs <- sum(!(x.ok <- is.finite(x)))) > 0) {
        is.not.finite.warning(x)
        x <- x[x.ok]
        warning(paste(bad.obs, "observations with NA/NaN/Inf in 'x' removed."))
    }
    x <- sort(x)
    n <- length(x)
    if (n < 2 || length(unique(x)) < 2) 
        stop("'x' must contain at least 2 non-missing distinct values.")
    pi.type <- match.arg(pi.type, c("upper", "lower"))
    if (pi.type == "upper") 
        lpl.rank <- 0
    else n.plus.one.minus.upl.rank <- 0
    if (!is.numeric(lpl.rank) || length(lpl.rank) > 1 || lpl.rank != 
        trunc(lpl.rank) || lpl.rank < 0 || lpl.rank >= n) 
        stop(paste("'lpl.rank' must be a non-negative integer less than", 
            "the number of non-missing finite values in 'x'"))
    if (pi.type == "lower" & lpl.rank < 1) 
        stop("When pi.type='lower', 'lpl.rank' must be a positive integer")
    if (!is.numeric(n.plus.one.minus.upl.rank) || length(n.plus.one.minus.upl.rank) > 
        1 || n.plus.one.minus.upl.rank != trunc(n.plus.one.minus.upl.rank) || 
        n.plus.one.minus.upl.rank < 0 || n.plus.one.minus.upl.rank >= 
        n) 
        stop(paste("'n.plus.one.minus.upl.rank' must be a non-negative integer", 
            "less than the number of non-missing finite values in 'x'"))
    if (pi.type == "upper" & n.plus.one.minus.upl.rank < 1) 
        stop("When pi.type='upper', 'n.plus.one.minus.upl.rank' must be a positive integer")
    upl.rank <- n + 1 - n.plus.one.minus.upl.rank
    args.list <- list(n = n, n.median = n.median, r = r, rule = rule, 
        lpl.rank = lpl.rank, n.plus.one.minus.upl.rank = n.plus.one.minus.upl.rank, 
        pi.type = pi.type, integrate.args.list = integrate.args.list)
    if (rule == "CA") 
        args.list <- c(args.list, list(m = m))
    if (rule == "k.of.m") 
        args.list <- c(args.list, list(k = k, m = m))
    conf.level <- do.call("predIntNparSimultaneousConfLevel", 
        args = args.list)
    upl.rank <- n + 1 - n.plus.one.minus.upl.rank
    switch(pi.type, lower = {
        lpl <- x[lpl.rank]
        upl <- ub
        lrl <- lpl.rank
        lru <- NULL
        n.plus.one.minus.upl.rank <- 0
    }, upper = {
        lpl <- lb
        upl <- x[upl.rank]
        lrl <- NULL
        lru <- upl.rank
        lpl.rank <- 0
    })
    ret.list <- list(distribution = "None", sample.size = n, 
        data.name = data.name, bad.obs = bad.obs)
    limits <- c(lpl, upl)
    names(limits) <- c("LPL", "UPL")
    string <- ifelse(rule == "k.of.m", "", paste("(", rule, " Rule)", 
        sep = ""))
    pi.obj <- list(name = "Prediction", rule = rule, limit.ranks = c(lrl, 
        lru), limits = limits, type = pi.type, method = paste("exact", 
        string), conf.level = conf.level, sample.size = n, k = k, 
        m = m, r = r, n.median = n.median)
    oldClass(pi.obj) <- "intervalEstimate"
    ret.list <- c(ret.list, list(interval = pi.obj))
    oldClass(ret.list) <- "estimate"
    ret.list
}
