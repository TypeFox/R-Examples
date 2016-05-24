predIntNparConfLevel <-
function (n, k = m, m = 1, lpl.rank = ifelse(pi.type == "upper", 
    0, 1), n.plus.one.minus.upl.rank = ifelse(pi.type == "lower", 
    0, 1), pi.type = "two.sided") 
{
    if (!is.vector(n, mode = "numeric") || !all(is.finite(n)) || 
        any(n != trunc(n)) || any(n < 2)) 
        stop("'n' must be a vector of positive integers greater than 1.")
    if (!is.vector(m, mode = "numeric") || !all(is.finite(m)) || 
        any(m != trunc(m)) || any(m < 1)) 
        stop("'m' must be a vector of positive integers")
    if (!is.vector(k, mode = "numeric") || !all(is.finite(k)) || 
        any(k != trunc(k)) || any(k < 1) || any(k > m)) 
        stop(paste("'k' must be a vector of positive integers,", 
            "and all values of 'k' must be between", "1 and the corresponding value of 'm'"))
    pi.type <- match.arg(pi.type, c("two.sided", "lower", "upper"))
    if (pi.type == "upper") 
        lpl.rank <- 0
    else if (pi.type == "lower") 
        n.plus.one.minus.upl.rank <- 0
    if (!is.vector(lpl.rank, mode = "numeric") || !all(is.finite(lpl.rank)) || 
        any(lpl.rank != trunc(lpl.rank)) || any(lpl.rank < 0 | 
        lpl.rank >= n)) 
        stop(paste("'lpl.rank' must be a vector of non-negative", 
            "integers less than the corresponding value of 'n'"))
    if (pi.type %in% c("two.sided", "lower") & any(lpl.rank < 
        1)) 
        stop("When pi.type='two.sided' or pi.type='lower', all values of 'lpl.rank' must be positive integers")
    if (!is.vector(n.plus.one.minus.upl.rank, mode = "numeric") || 
        !all(is.finite(n.plus.one.minus.upl.rank)) || any(n.plus.one.minus.upl.rank != 
        trunc(n.plus.one.minus.upl.rank)) || any(n.plus.one.minus.upl.rank < 
        0 | n.plus.one.minus.upl.rank >= n)) 
        stop(paste("'n.plus.one.minus.upl.rank' must be a vector of non-negative", 
            "integers less than the corresponding value of 'n'"))
    if (pi.type %in% c("two.sided", "upper") & any(n.plus.one.minus.upl.rank < 
        1)) 
        stop("When pi.type='two.sided' or pi.type='upper' all values of 'n.plus.one.minus.upl.rank' must be positive integers")
    arg.mat <- cbind.no.warn(n = as.vector(n), k = as.vector(k), 
        m = as.vector(m), lpl.rank = as.vector(lpl.rank), n.plus.one.minus.upl.rank = as.vector(n.plus.one.minus.upl.rank))
    for (i in c("n", "k", "m", "lpl.rank", "n.plus.one.minus.upl.rank")) assign(i, 
        arg.mat[, i])
    if (pi.type == "two.sided") {
        upl.rank <- n + 1 - n.plus.one.minus.upl.rank
        if (any(lpl.rank >= upl.rank)) 
            stop(paste("Illegal values for 'lpl.rank' and 'n.plus.one.minus.upl.rank'.", 
                "Make one or both of them smaller"))
    }
    N <- length(n)
    conf.level.vec <- numeric(N)
    for (i in 1:N) {
        conf.level.vec[i] <- predIntNparConfLevelScalar(n = n[i], 
            k = k[i], m = m[i], lpl.rank = lpl.rank[i], n.plus.one.minus.upl.rank = n.plus.one.minus.upl.rank[i], 
            pi.type = pi.type)
    }
    names(conf.level.vec) <- NULL
    conf.level.vec
}
