ciNparConfLevelScalar <-
function (n, p = 0.5, lcl.rank = ifelse(ci.type == "upper", 0, 
    1), n.plus.one.minus.ucl.rank = ifelse(ci.type == "lower", 
    0, 1), ci.type = c("two.sided", "lower", "upper")) 
{
    ci.type <- match.arg(ci.type)
    if (!is.numeric(n) || length(n) > 1 || n != trunc(n) || n < 
        2) 
        stop("'n' must be a positive integer greater than 1.")
    if (!is.numeric(p) || length(p) > 1 || p <= 0 || p >= 1) 
        stop("'p' must be greater than 0 and less than 1.")
    if (!is.numeric(lcl.rank) || length(lcl.rank) > 1 || lcl.rank != 
        trunc(lcl.rank) || lcl.rank < 0) 
        stop("'lcl.rank' must be a non-negative integer")
    if (!is.numeric(n.plus.one.minus.ucl.rank) || length(n.plus.one.minus.ucl.rank) > 
        1 || n.plus.one.minus.ucl.rank != trunc(n.plus.one.minus.ucl.rank) || 
        n.plus.one.minus.ucl.rank < 0) 
        stop("'n.plus.one.minus.ucl.rank' must be a non-negative integer")
    ucl.rank <- (n + 1) - n.plus.one.minus.ucl.rank
    if (ci.type == "two.sided" && lcl.rank >= ucl.rank) 
        stop(paste("Illegal values for 'lcl.rank' and 'n.plus.one.minus.ucl.rank'. ", 
            "Make one or both of them smaller"))
    switch(ci.type, two.sided = {
        if (lcl.rank == 0 || n.plus.one.minus.ucl.rank == 0) stop(paste("'lcl.rank' and 'n.plus.one.minus.ucl.rank'", 
            "must both be at least 1 when ci.type='two.sided'"))
    }, lower = {
        if (lcl.rank == 0) stop(paste("'lcl.rank' must be at least 1", 
            "when ci.type='lower'"))
        n.plus.one.minus.ucl.rank <- 0
    }, upper = {
        if (n.plus.one.minus.ucl.rank == 0) stop(paste("'n.plus.one.minus.ucl.rank' must be", 
            "at least 1 when ci.type='upper'"))
        lcl.rank <- 0
    })
    conf.level <- switch(ci.type, two.sided = pbinom(ucl.rank - 
        1, n, p) - pbinom(lcl.rank - 1, n, p), lower = 1 - pbinom(lcl.rank - 
        1, n, p), upper = pbinom(ucl.rank - 1, n, p))
    names(conf.level) <- NULL
    conf.level
}
