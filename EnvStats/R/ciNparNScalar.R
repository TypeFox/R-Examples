ciNparNScalar <-
function (p = 0.5, lcl.rank = ifelse(ci.type == "upper", 0, 1), 
    n.plus.one.minus.ucl.rank = ifelse(ci.type == "lower", 0, 
        1), ci.type = c("two.sided", "lower", "upper"), conf.level = 0.95) 
{
    if (!is.numeric(p) || length(p) > 1 || p <= 0 || p >= 1) 
        stop("'p' must be a scalar greater than 0 and less than 1.")
    ci.type <- match.arg(ci.type)
    if (!is.numeric(lcl.rank) || length(lcl.rank) > 1 || lcl.rank != 
        trunc(lcl.rank) || lcl.rank < 0) 
        stop("'lcl.rank' must be a non-negative integer")
    if (!is.numeric(n.plus.one.minus.ucl.rank) || length(n.plus.one.minus.ucl.rank) > 
        1 || n.plus.one.minus.ucl.rank != trunc(n.plus.one.minus.ucl.rank) || 
        n.plus.one.minus.ucl.rank < 0) 
        stop("'n.plus.one.minus.ucl.rank' must be a non-negative integer")
    if (!is.numeric(conf.level) || length(conf.level) > 1 || 
        conf.level <= 0 || conf.level >= 1) 
        stop("'conf.level' must be a scalar greater than 0 and less than 1.")
    if (ci.type == "lower") {
        fcn.to.min.lower <- function(n, p.weird, lcl.rank, conf.level) {
            ((1 - pbeta(1 - p.weird, n - (lcl.rank - 1), lcl.rank)) - 
                conf.level)^2
        }
        start <- 2
        n <- ceiling(nlminb(start = start, objective = fcn.to.min.lower, 
            lower = 2, p.weird = p, lcl.rank = lcl.rank, conf.level = conf.level)$par)
    }
    else if (ci.type == "upper") {
        fcn.to.min.upper <- function(n, p.weird, n.plus.one.minus.ucl.rank, 
            conf.level) {
            ucl.rank <- n + 1 - n.plus.one.minus.ucl.rank
            (pbeta(1 - p.weird, n - (ucl.rank - 1), ucl.rank) - 
                conf.level)^2
        }
        start <- 2
        n <- ceiling(nlminb(start = start, objective = fcn.to.min.upper, 
            lower = 2, p.weird = p, n.plus.one.minus.ucl.rank = n.plus.one.minus.ucl.rank, 
            conf.level = conf.level)$par)
    }
    else {
        fcn.to.min <- function(n, p.weird, lcl.rank, n.plus.one.minus.ucl.rank, 
            conf.level) {
            ucl.rank <- n + 1 - n.plus.one.minus.ucl.rank
            ((pbeta(1 - p.weird, n - (ucl.rank - 1), ucl.rank) - 
                pbeta(1 - p.weird, n - (lcl.rank - 1), lcl.rank)) - 
                conf.level)^2
        }
        start <- max(2, lcl.rank + n.plus.one.minus.ucl.rank)
        n <- ceiling(nlminb(start = start, objective = fcn.to.min, 
            lower = 2, p.weird = p, lcl.rank = lcl.rank, n.plus.one.minus.ucl.rank = n.plus.one.minus.ucl.rank, 
            conf.level = conf.level)$par)
    }
    if (n > start) {
        test.conf.level <- ciNparConfLevelScalar(n = n - 1, p = p, 
            lcl.rank = lcl.rank, n.plus.one.minus.ucl.rank = n.plus.one.minus.ucl.rank, 
            ci.type = ci.type)
        if (test.conf.level >= conf.level) 
            n <- n - 1
    }
    n
}
