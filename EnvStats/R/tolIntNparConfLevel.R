tolIntNparConfLevel <-
function (n, coverage = 0.95, ltl.rank = ifelse(ti.type == "upper", 
    0, 1), n.plus.one.minus.utl.rank = ifelse(ti.type == "lower", 
    0, 1), ti.type = "two.sided") 
{
    if (!is.vector(n, mode = "numeric") || !all(is.finite(n)) || 
        any(n != trunc(n)) || any(n < 2)) 
        stop("'n' must be a vector of positive integers greater than 1.")
    if (!is.vector(coverage, mode = "numeric") || !all(is.finite(coverage)) || 
        any(coverage <= 0) || any(coverage >= 1)) 
        stop("All values of 'coverage' must be greater than 0 and less than 1.")
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
    arg.mat <- cbind.no.warn(n = as.vector(n), coverage = as.vector(coverage), 
        ltl.rank = as.vector(ltl.rank), n.plus.one.minus.utl.rank = as.vector(n.plus.one.minus.utl.rank))
    for (i in c("n", "coverage", "ltl.rank", "n.plus.one.minus.utl.rank")) assign(i, 
        arg.mat[, i])
    if (ti.type == "two.sided") {
        utl.rank <- n + 1 - n.plus.one.minus.utl.rank
        if (any(ltl.rank >= utl.rank)) 
            stop(paste("Illegal values for 'ltl.rank' and 'n.plus.one.minus.utl.rank'.", 
                "Make one or both of them smaller"))
    }
    N <- length(n)
    conf.level.vec <- numeric(N)
    for (i in 1:N) {
        conf.level.vec[i] <- tolIntNparConfLevelScalar(n = n[i], 
            coverage = coverage[i], ltl.rank = ltl.rank[i], n.plus.one.minus.utl.rank = n.plus.one.minus.utl.rank[i], 
            ti.type = ti.type)
    }
    names(conf.level.vec) <- NULL
    conf.level.vec
}
