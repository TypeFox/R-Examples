tolIntNpar <-
function (x, coverage, conf.level, cov.type = "content", ltl.rank = ifelse(ti.type == 
    "upper", 0, 1), n.plus.one.minus.utl.rank = ifelse(ti.type == 
    "lower", 0, 1), lb = -Inf, ub = Inf, ti.type = "two-sided") 
{
    if (!is.vector(x, mode = "numeric")) 
        stop("'x' must be a numeric vector.")
    cov.type <- match.arg(cov.type, c("content", "expectation"))
    if (cov.type == "content") {
        miss.cov <- missing(coverage)
        miss.con <- missing(conf.level)
        if ((miss.cov && miss.con) || (!miss.cov && !miss.con)) 
            stop(paste("You must supply a value for 'coverage'", 
                "or 'conf.level', but not both"))
        if (!miss.cov) {
            if (!is.numeric(coverage) || length(coverage) > 1 || 
                coverage <= 0 || coverage >= 1) 
                stop("'coverage' must be a scalar greater than 0 and less than 1.")
        }
        else {
            if (!is.numeric(conf.level) || length(conf.level) > 
                1 || conf.level <= 0 || conf.level >= 1) 
                stop("'conf.level' must be a scalar greater than 0 and less than 1.")
        }
    }
    ti.type <- match.arg(ti.type, c("two-sided", "lower", "upper"))
    if (!is.numeric(ltl.rank) || length(ltl.rank) > 1 || ltl.rank != 
        trunc(ltl.rank) || ltl.rank < 0 || !is.numeric(n.plus.one.minus.utl.rank) || 
        length(n.plus.one.minus.utl.rank) > 1 || n.plus.one.minus.utl.rank != 
        trunc(n.plus.one.minus.utl.rank) || n.plus.one.minus.utl.rank < 
        0) 
        stop("'ltl.rank' and 'n.plus.one.minus.utl.rank' must be non-negative integers")
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
    utl.rank <- n + 1 - n.plus.one.minus.utl.rank
    if (ltl.rank >= utl.rank) 
        stop(paste("Illegal values for 'ltl.rank' and 'n.plus.one.minus.utl.rank'. ", 
            "Make one or both of them smaller"))
    switch(ti.type, `two-sided` = {
        if (ltl.rank == 0 || n.plus.one.minus.utl.rank == 0) stop(paste("'ltl.rank' and 'n.plus.one.minus.utl.rank'", 
            "must both be at least 1 when ti.type='two-sided'"))
        ltl <- x[ltl.rank]
        utl <- x[utl.rank]
        lrl <- ltl.rank
        lru <- utl.rank
    }, lower = {
        if (ltl.rank == 0) stop(paste("'ltl.rank' must be at least 1", 
            "when ti.type='lower'"))
        ltl <- x[ltl.rank]
        utl <- ub
        lrl <- ltl.rank
        lru <- NULL
        n.plus.one.minus.utl.rank <- 0
    }, upper = {
        if (n.plus.one.minus.utl.rank == 0) stop(paste("'n.plus.one.minus.utl.rank' must be", 
            "at least 1 when ti.type='upper'"))
        ltl <- lb
        utl <- x[utl.rank]
        lrl <- NULL
        lru <- utl.rank
        ltl.rank <- 0
    })
    ret.list <- list(distribution = "None", sample.size = n, 
        data.name = data.name, bad.obs = bad.obs)
    limits <- c(ltl, utl)
    names(limits) <- c("LTL", "UTL")
    if (cov.type == "content") {
        if (!miss.cov) 
            conf.level <- 1 - pbeta(coverage, utl.rank - ltl.rank, 
                n.plus.one.minus.utl.rank + ltl.rank)
        else coverage <- qbeta(1 - conf.level, utl.rank - ltl.rank, 
            n.plus.one.minus.utl.rank + ltl.rank)
        ti.obj <- list(name = "Tolerance", coverage = coverage, 
            coverage.type = cov.type, limit.ranks = c(lrl, lru), 
            limits = limits, type = ti.type, method = "Exact", 
            conf.level = conf.level, sample.size = n)
    }
    else {
        coverage <- (utl.rank - ltl.rank)/(n + 1)
        ti.obj <- list(name = "Tolerance", coverage = coverage, 
            coverage.type = cov.type, limit.ranks = c(lrl, lru), 
            limits = limits, type = ti.type, method = "Exact", 
            sample.size = n)
    }
    oldClass(ti.obj) <- "intervalEstimate"
    ret.list <- c(ret.list, list(interval = ti.obj))
    oldClass(ret.list) <- "estimate"
    ret.list
}
