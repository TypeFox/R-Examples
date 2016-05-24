ciNparConfLevel <-
function (n, p = 0.5, lcl.rank = ifelse(ci.type == "upper", 0, 
    1), n.plus.one.minus.ucl.rank = ifelse(ci.type == "lower", 
    0, 1), ci.type = "two.sided") 
{
    ci.type <- match.arg(ci.type, c("two.sided", "lower", "upper"))
    if (!is.vector(n, mode = "numeric") || !all(is.finite(n)) || 
        any(n != trunc(n)) || any(n < 2)) 
        stop("'n' must be a vector of positive integers greater than 1.")
    if (!is.vector(p, mode = "numeric") || !all(is.finite(p)) || 
        any(p <= 0) || any(p >= 1)) 
        stop("All values of 'p' must be greater than 0 and less than 1.")
    if (!is.vector(lcl.rank, mode = "numeric") || !all(is.finite(lcl.rank)) || 
        any(lcl.rank != trunc(lcl.rank)) || any(lcl.rank < 0)) 
        stop("'lcl.rank' must be a vector of nonnegative integers")
    if (!is.vector(n.plus.one.minus.ucl.rank, mode = "numeric") || 
        !all(is.finite(n.plus.one.minus.ucl.rank)) || any(n.plus.one.minus.ucl.rank != 
        trunc(n.plus.one.minus.ucl.rank)) || any(n.plus.one.minus.ucl.rank < 
        0)) 
        stop("'n.plus.one.minus.ucl.rank' must be a vector of nonnegative integers")
    arg.mat <- cbind.no.warn(n = as.vector(n), p = as.vector(p), 
        lcl.rank = as.vector(lcl.rank), n.plus.one.minus.ucl.rank = as.vector(n.plus.one.minus.ucl.rank))
    for (i in c("n", "p", "lcl.rank", "n.plus.one.minus.ucl.rank")) assign(i, 
        arg.mat[, i])
    N <- length(n)
    conf.level.vec <- numeric(N)
    for (i in 1:N) {
        conf.level.vec[i] <- ciNparConfLevelScalar(n = n[i], 
            p = p[i], lcl.rank = lcl.rank[i], n.plus.one.minus.ucl.rank = n.plus.one.minus.ucl.rank[i], 
            ci.type = ci.type)
    }
    names(conf.level.vec) <- NULL
    conf.level.vec
}
