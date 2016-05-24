tolIntNparConfLevelScalar <-
function (n, coverage = 0.95, ltl.rank = ifelse(ti.type == "upper", 
    0, 1), n.plus.one.minus.utl.rank = ifelse(ti.type == "lower", 
    0, 1), ti.type = "two.sided") 
{
    if (!is.numeric(n) || length(n) > 1 || n != trunc(n) || n < 
        2) 
        stop("'n' must be an integer greater than 1.")
    if (!is.numeric(coverage) || length(coverage) > 1 || coverage <= 
        0 || coverage >= 1) 
        stop("'coverage' must be greater than 0 and less than 1.")
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
    conf.level <- 1 - pbeta(coverage, utl.rank - ltl.rank, n.plus.one.minus.utl.rank + 
        ltl.rank)
    names(conf.level) <- NULL
    conf.level
}
