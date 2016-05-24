predIntNparNScalar <-
function (k = m, m = 1, lpl.rank = ifelse(pi.type == "upper", 
    0, 1), n.plus.one.minus.upl.rank = ifelse(pi.type == "lower", 
    0, 1), pi.type = c("two.sided", "lower", "upper"), conf.level = 0.95, 
    n.max = 5000, maxiter = 1000) 
{
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
        trunc(lpl.rank) || lpl.rank < 0) 
        stop("'lpl.rank' must be a non-negative integer")
    if (pi.type %in% c("two.sided", "lower") & lpl.rank < 1) 
        stop("When pi.type='two.sided' or pi.type='lower', 'lpl.rank' must be a positive integer")
    if (!is.numeric(n.plus.one.minus.upl.rank) || length(n.plus.one.minus.upl.rank) > 
        1 || n.plus.one.minus.upl.rank != trunc(n.plus.one.minus.upl.rank) || 
        n.plus.one.minus.upl.rank < 0) 
        stop("'n.plus.one.minus.upl.rank' must be a non-negative integer")
    if (pi.type %in% c("two.sided", "upper") & n.plus.one.minus.upl.rank < 
        1) 
        stop("When pi.type='two.sided' or pi.type='upper', 'n.plus.one.minus.upl.rank' must be a positive integer")
    n.min <- max(2, lpl.rank + n.plus.one.minus.upl.rank)
    if (n.max < n.min) 
        stop("Argument 'n.max' is too small relative to 'lpl.rank' and 'n.plus.one.minus.upl.rank'")
    conf.level.n.min <- predIntNparConfLevelScalar(n = n.min, 
        k = k, m = m, lpl.rank = lpl.rank, n.plus.one.minus.upl.rank = n.plus.one.minus.upl.rank, 
        pi.type = pi.type)
    if (conf.level.n.min >= conf.level) {
        n <- n.min
    }
    else {
        conf.level.n.max <- predIntNparConfLevelScalar(n = n.max, 
            k = k, m = m, lpl.rank = lpl.rank, n.plus.one.minus.upl.rank = n.plus.one.minus.upl.rank, 
            pi.type = pi.type)
        if (conf.level.n.max < conf.level) {
            n <- NA
            warning("Error in algorithm.  Try increasing the value of 'n.max'")
        }
        else {
            simple <- k == m
            if (pi.type == "two.sided") {
                a <- 1 - conf.level
                b <- conf.level * (1 - 2 * m) - 1
                cc <- conf.level * m * (1 - m)
                n <- trunc((-b + sqrt(b^2 - 4 * a * cc))/(2 * 
                  a))
                if ((n * (n - 1))/((n + m) * (n + m - 1)) < conf.level) 
                  n <- n + 1
                if (lpl.rank > 1 || n.plus.one.minus.upl.rank > 
                  1) 
                  simple <- FALSE
            }
            else {
                n <- trunc((m * conf.level)/(1 - conf.level))
                if (n/(n + m) < conf.level) 
                  n <- n + 1
                if ((pi.type == "upper" && n.plus.one.minus.upl.rank > 
                  1) || (pi.type == "lower" && lpl.rank > 1)) {
                  simple <- FALSE
                }
            }
            if (!simple) {
                fcn.for.root <- function(n.weird, m.weird, k.weird, 
                  lpl.rank.weird, n.plus.one.minus.upl.rank.weird, 
                  conf.level) {
                  n <- trunc(n.weird)
                  m <- m.weird
                  k <- k.weird
                  lpl.rank <- lpl.rank.weird
                  n.plus.one.minus.upl.rank <- n.plus.one.minus.upl.rank.weird
                  vec <- k:m
                  sum(exp(lchoose(m - vec + lpl.rank + n.plus.one.minus.upl.rank - 
                    1, m - vec) + lchoose(vec + n - lpl.rank - 
                    n.plus.one.minus.upl.rank, vec) - lchoose(n + 
                    m, m))) - conf.level
                }
                uni.list <- uniroot(fcn.for.root, interval = c(n.min, 
                  n.max), m.weird = m, k.weird = k, lpl.rank.weird = lpl.rank, 
                  n.plus.one.minus.upl.rank.weird = n.plus.one.minus.upl.rank, 
                  conf.level = conf.level, tol = 1, maxiter = maxiter)
                if (uni.list$f.root > 0) 
                  n <- trunc(uni.list$root)
                else n <- trunc(uni.list$root) + 1
            }
        }
    }
    n
}
