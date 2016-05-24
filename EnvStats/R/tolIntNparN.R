tolIntNparN <-
function (coverage = 0.95, conf.level = 0.95, cov.type = "content", 
    ltl.rank = ifelse(ti.type == "upper", 0, 1), n.plus.one.minus.utl.rank = ifelse(ti.type == 
        "lower", 0, 1), ti.type = "two.sided") 
{
    ti.type <- match.arg(ti.type, c("two.sided", "lower", "upper"))
    if (ti.type == "upper") 
        ltl.rank <- 0
    else if (ti.type == "lower") 
        n.plus.one.minus.utl.rank <- 0
    if (!is.vector(ltl.rank, mode = "numeric") || !all(is.finite(ltl.rank)) || 
        any(ltl.rank != trunc(ltl.rank)) || any(ltl.rank < 0)) 
        stop("'ltl.rank' must be a vector of non-negative integers")
    if (ti.type %in% c("two.sided", "lower") & any(ltl.rank < 
        1)) 
        stop("When ti.type='two.sided' or ti.type='lower', all values of 'ltl.rank' must be positive integers")
    if (!is.vector(n.plus.one.minus.utl.rank, mode = "numeric") || 
        !all(is.finite(n.plus.one.minus.utl.rank)) || any(n.plus.one.minus.utl.rank != 
        trunc(n.plus.one.minus.utl.rank)) || any(n.plus.one.minus.utl.rank < 
        0)) 
        stop("'n.plus.one.minus.utl.rank' must be a vector of non-negative integers")
    if (ti.type %in% c("two.sided", "upper") & any(n.plus.one.minus.utl.rank < 
        1)) 
        stop("When ti.type='two.sided' or ti.type='upper' all values of 'n.plus.one.minus.utl.rank' must be positive integers")
    if (!all(is.finite(coverage)) || any(coverage <= 0) || any(coverage >= 
        1)) 
        stop("All values of 'coverage' must be greater than 0 and less than 1.")
    cov.type <- match.arg(cov.type, c("content", "expectation"))
    if (!all(is.finite(conf.level)) || any(conf.level <= 0) || 
        any(conf.level >= 1)) 
        stop("All values of 'conf.level' must be greater than 0 and less than 1.")
    arg.mat <- cbind.no.warn(ltl.rank = as.vector(ltl.rank), 
        n.plus.one.minus.utl.rank = as.vector(n.plus.one.minus.utl.rank), 
        coverage = as.vector(coverage), conf.level = as.vector(conf.level))
    for (i in c("ltl.rank", "n.plus.one.minus.utl.rank", "coverage", 
        "conf.level")) assign(i, arg.mat[, i])
    N <- length(ltl.rank)
    n.vec <- numeric(N)
    for (i in 1:N) {
        n.vec[i] <- tolIntNparNScalar(ltl.rank = ltl.rank[i], 
            n.plus.one.minus.utl.rank = n.plus.one.minus.utl.rank[i], 
            coverage = coverage[i], cov.type = cov.type, ti.type = ti.type, 
            conf.level = conf.level[i])
    }
    names(n.vec) <- NULL
    n.vec
}
