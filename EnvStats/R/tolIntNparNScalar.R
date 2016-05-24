tolIntNparNScalar <-
function (ltl.rank = ifelse(ti.type == "upper", 0, 1), n.plus.one.minus.utl.rank = ifelse(ti.type == 
    "lower", 0, 1), coverage = 0.95, cov.type = "content", ti.type = "two.sided", 
    conf.level = 0.95) 
{
    ti.type <- match.arg(ti.type, c("two.sided", "lower", "upper"))
    if (ti.type == "upper" & ltl.rank != 0) 
        ltl.rank <- 0
    else if (ti.type == "lower" & n.plus.one.minus.utl.rank != 
        0) 
        n.plus.one.minus.utl.rank <- 0
    if (!is.numeric(ltl.rank) || length(ltl.rank) > 1 || ltl.rank != 
        trunc(ltl.rank)) 
        stop("'ltl.rank' must be an integer")
    if (ti.type == "two.sided" & ltl.rank < 1) 
        stop("When ti.type='two.sided', 'ltl.rank' must be a positive integer")
    if (!is.numeric(n.plus.one.minus.utl.rank) || length(n.plus.one.minus.utl.rank) > 
        1 || n.plus.one.minus.utl.rank != trunc(n.plus.one.minus.utl.rank)) 
        stop("'n.plus.one.minus.utl.rank' must be an integer")
    if (ti.type == "two.sided" & n.plus.one.minus.utl.rank < 
        1) 
        stop("When ti.type='two.sided', 'n.plus.one.minus.utl.rank' must be a positive integer")
    if (!is.numeric(coverage) || length(coverage) > 1 || coverage <= 
        0 || coverage >= 1) 
        stop("'coverage' must be a scalar greater than 0 and less than 1.")
    cov.type <- match.arg(cov.type, c("content", "expectation"))
    if (!is.numeric(conf.level) || length(conf.level) > 1 || 
        conf.level <= 0 || conf.level >= 1) 
        stop("'conf.level' must be a scalar greater than 0 and less than 1.")
    alpha <- 1 - conf.level
    if (cov.type == "content") {
        start <- round((0.25 * qchisq(conf.level, 2 * (ltl.rank + 
            n.plus.one.minus.utl.rank)) * (1 + coverage))/(1 - 
            coverage) + 0.5 * (ltl.rank + n.plus.one.minus.utl.rank - 
            1), 0)
        fcn.to.min <- function(n, ltl.rank, n.plus.one.minus.utl.rank, 
            alpha, coverage) {
            (pbeta(coverage, n + 1 - n.plus.one.minus.utl.rank - 
                ltl.rank, n.plus.one.minus.utl.rank + ltl.rank) - 
                alpha)^2
        }
        n <- ceiling(nlminb(start = start, objective = fcn.to.min, 
            lower = 2, ltl.rank = ltl.rank, n.plus.one.minus.utl.rank = n.plus.one.minus.utl.rank, 
            alpha = alpha, coverage = coverage)$par)
    }
    else {
        n <- ceiling((ltl.rank + n.plus.one.minus.utl.rank)/(1 - 
            coverage) - 1)
    }
    n
}
