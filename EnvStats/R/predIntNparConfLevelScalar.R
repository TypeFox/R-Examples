predIntNparConfLevelScalar <-
function (n, k = m, m = 1, lpl.rank = ifelse(pi.type == "upper", 
    0, 1), n.plus.one.minus.upl.rank = ifelse(pi.type == "lower", 
    0, 1), pi.type = c("two.sided", "lower", "upper")) 
{
    if (!is.numeric(n) || length(n) > 1 || n != trunc(n) || n < 
        2) 
        stop("'n' must be a positive integer greater than 1.")
    if (!is.numeric(m) || length(m) > 1 || m != trunc(m) || m < 
        1) 
        stop("'m' must be a positive integer")
    if (!is.numeric(k) || length(k) > 1 || k != trunc(k) || k < 
        1 || k > m) 
        stop("'k' must be a positive integer between 1 and 'm'")
    pi.type <- match.arg(pi.type)
    if (pi.type == "upper") 
        lpl.rank <- 0
    else if (pi.type == "lower") 
        n.plus.one.minus.upl.rank <- 0
    if (!is.numeric(lpl.rank) || length(lpl.rank) > 1 || lpl.rank != 
        trunc(lpl.rank) || lpl.rank < 0 || lpl.rank >= n) 
        stop("'lpl.rank' must be a non-negative integer less than 'n'")
    if (pi.type %in% c("two.sided", "lower") & lpl.rank < 1) 
        stop(paste("When pi.type='two.sided' or pi.type='lower',", 
            "'lpl.rank' must be a positive integer"))
    if (!is.numeric(n.plus.one.minus.upl.rank) || length(n.plus.one.minus.upl.rank) > 
        1 || n.plus.one.minus.upl.rank != trunc(n.plus.one.minus.upl.rank) || 
        n.plus.one.minus.upl.rank < 0 || n.plus.one.minus.upl.rank >= 
        n) 
        stop("'n.plus.one.minus.upl.rank' must be a non-negative integer less than 'n'")
    if (pi.type %in% c("two.sided", "upper") & n.plus.one.minus.upl.rank < 
        1) 
        stop(paste("When pi.type='two.sided' or pi.type='upper',", 
            "'n.plus.one.minus.upl.rank' must be a positive integer"))
    upl.rank <- n + 1 - n.plus.one.minus.upl.rank
    if (pi.type == "two.sided" && lpl.rank >= upl.rank) 
        stop(paste("Illegal values for 'lpl.rank' and 'n.plus.one.minus.upl.rank'. ", 
            "Make one or both of them smaller."))
    vec <- k:m
    conf.level <- sum(exp(lchoose(m - vec + lpl.rank + n.plus.one.minus.upl.rank - 
        1, m - vec) + lchoose(vec + n - lpl.rank - n.plus.one.minus.upl.rank, 
        vec) - lchoose(n + m, m)))
    names(conf.level) <- NULL
    conf.level
}
