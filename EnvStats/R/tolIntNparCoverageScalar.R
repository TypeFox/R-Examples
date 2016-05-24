tolIntNparCoverageScalar <-
function (n, conf.level = 0.95, cov.type = "content", ltl.rank = ifelse(ti.type == 
    "upper", 0, 1), n.plus.one.minus.utl.rank = ifelse(ti.type == 
    "lower", 0, 1), ti.type = "two.sided") 
{
    if (!is.numeric(n) || length(n) > 1 || n != trunc(n) || n < 
        2) 
        stop("'n' must be a positive integer greater than 1.")
    cov.type <- match.arg(cov.type, c("content", "expectation"))
    if (cov.type == "content") {
        if (!is.numeric(conf.level) || length(conf.level) > 1 || 
            conf.level <= 0 || conf.level >= 1) 
            stop("'conf.level' must be a scalar greater than 0 and less than 1.")
    }
    ti.type <- match.arg(ti.type, c("two.sided", "lower", "upper"))
    if (!is.numeric(ltl.rank) || length(ltl.rank) > 1 || ltl.rank != 
        trunc(ltl.rank) || ltl.rank < 0) 
        stop("'ltl.rank' must be a non-negative integer")
    if (!is.numeric(n.plus.one.minus.utl.rank) || length(n.plus.one.minus.utl.rank) > 
        1 || n.plus.one.minus.utl.rank != trunc(n.plus.one.minus.utl.rank) || 
        n.plus.one.minus.utl.rank < 0) 
        stop("'n.plus.one.minus.utl.rank' must be a non-negative integer")
    utl.rank <- n + 1 - n.plus.one.minus.utl.rank
    if (ltl.rank >= utl.rank) 
        stop(paste("Illegal values for 'ltl.rank' and 'n.plus.one.minus.utl.rank'. ", 
            "Make one or both of them smaller"))
    switch(ti.type, two.sided = {
        if (ltl.rank == 0 || n.plus.one.minus.utl.rank == 0) stop(paste("'ltl.rank' and 'n.plus.one.minus.utl.rank'", 
            "must both be at least 1 when ti.type='two.sided'"))
    }, lower = {
        if (ltl.rank == 0) stop(paste("'ltl.rank' must be at least 1", 
            "when ti.type='lower'"))
        n.plus.one.minus.utl.rank <- 0
    }, upper = {
        if (n.plus.one.minus.utl.rank == 0) stop(paste("'n.plus.one.minus.utl.rank' must be", 
            "at least 1 when ti.type='upper'"))
        ltl.rank <- 0
    })
    if (cov.type == "content") {
        coverage <- qbeta(1 - conf.level, utl.rank - ltl.rank, 
            n.plus.one.minus.utl.rank + ltl.rank)
    }
    else {
        coverage <- (utl.rank - ltl.rank)/(n + 1)
    }
    names(coverage) <- NULL
    coverage
}
