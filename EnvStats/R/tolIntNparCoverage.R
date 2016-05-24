tolIntNparCoverage <-
function (n, conf.level = 0.95, cov.type = "content", ltl.rank = ifelse(ti.type == 
    "upper", 0, 1), n.plus.one.minus.utl.rank = ifelse(ti.type == 
    "lower", 0, 1), ti.type = "two.sided") 
{
    if (!is.vector(n, mode = "numeric") || !all(is.finite(n)) || 
        any(n != trunc(n)) || any(n < 2)) 
        stop("'n' must be a vector of positive integers greater than 1.")
    cov.type <- match.arg(cov.type, c("content", "expectation"))
    if (cov.type == "content") {
        if (!is.vector(conf.level, mode = "numeric") || !all(is.finite(conf.level)) || 
            any(conf.level <= 0) || any(conf.level >= 1)) 
            stop("All values of 'conf.level' must be greater than 0 and less than 1.")
    }
    ti.type <- match.arg(ti.type, c("two.sided", "lower", "upper"))
    if (ti.type == "upper") 
        ltl.rank <- 0
    else if (ti.type == "lower") 
        n.plus.one.minus.utl.rank <- 0
    if (!is.vector(ltl.rank, mode = "numeric") || !all(is.finite(ltl.rank)) || 
        any(ltl.rank != trunc(ltl.rank)) || any(ltl.rank < 0 | 
        ltl.rank >= n)) 
        stop(paste("'ltl.rank' must be a vector of non-negative", 
            "integers less than the corresponding value of 'n'"))
    if (ti.type %in% c("two.sided", "lower") & any(ltl.rank < 
        1)) 
        stop("When ti.type='two.sided' or ti.type='lower', all values of 'ltl.rank' must be positive integers")
    if (!is.vector(n.plus.one.minus.utl.rank, mode = "numeric") || 
        !all(is.finite(n.plus.one.minus.utl.rank)) || any(n.plus.one.minus.utl.rank != 
        trunc(n.plus.one.minus.utl.rank)) || any(n.plus.one.minus.utl.rank < 
        0 | n.plus.one.minus.utl.rank >= n)) 
        stop(paste("'n.plus.one.minus.utl.rank' must be a vector of non-negative", 
            "integers less than the corresponding value of 'n'"))
    if (ti.type %in% c("two.sided", "upper") & any(n.plus.one.minus.utl.rank < 
        1)) 
        stop("When ti.type='two.sided' or ti.type='upper' all values of 'n.plus.one.minus.utl.rank' must be positive integers")
    arg.mat <- cbind.no.warn(n = as.vector(n), conf.level = as.vector(conf.level), 
        ltl.rank = as.vector(ltl.rank), n.plus.one.minus.utl.rank = as.vector(n.plus.one.minus.utl.rank))
    for (i in c("n", "conf.level", "ltl.rank", "n.plus.one.minus.utl.rank")) assign(i, 
        arg.mat[, i])
    if (ti.type == "two.sided") {
        utl.rank <- n + 1 - n.plus.one.minus.utl.rank
        if (any(ltl.rank >= utl.rank)) 
            stop(paste("Illegal values for 'ltl.rank' and 'n.plus.one.minus.utl.rank'.", 
                "Make one or both of them smaller"))
    }
    N <- length(n)
    coverage.vec <- numeric(N)
    for (i in 1:N) {
        coverage.vec[i] <- tolIntNparCoverageScalar(n = n[i], 
            conf.level = conf.level[i], cov.type = cov.type, 
            ltl.rank = ltl.rank[i], n.plus.one.minus.utl.rank = n.plus.one.minus.utl.rank[i], 
            ti.type = ti.type)
    }
    names(coverage.vec) <- NULL
    coverage.vec
}
