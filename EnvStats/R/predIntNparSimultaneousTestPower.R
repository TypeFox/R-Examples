predIntNparSimultaneousTestPower <-
function (n, n.median = 1, k = 1, m = 2, r = 1, rule = "k.of.m", 
    lpl.rank = ifelse(pi.type == "upper", 0, 1), n.plus.one.minus.upl.rank = ifelse(pi.type == 
        "lower", 0, 1), delta.over.sigma = 0, pi.type = "upper", 
    r.shifted = r, method = "approx", NMC = 100, ci = FALSE, 
    ci.conf.level = 0.95, integrate.args.list = NULL) 
{
    rule <- match.arg(rule, c("k.of.m", "CA", "Modified.CA"), 
        several.ok = TRUE)
    pi.type <- match.arg(pi.type, c("upper", "lower"))
    method <- match.arg(method, c("approx", "simulate"))
    if (!is.vector(n, mode = "numeric") || !all(is.finite(n)) || 
        any(n < 2)) 
        stop(paste("'n' must be a numeric vector", "with all elements greater than or equal to 2", 
            "and no Missing (NA), Infinite (-Inf, Inf),", "or Undefined (Nan) values."))
    if (!is.vector(n.median, mode = "numeric") || !all(is.finite(n.median)) || 
        any(n.median < 1) || !all(n.median == trunc(n.median)) || 
        !all(is.odd(n.median))) 
        stop("'n.median' must be a numeric vector of positive odd integers")
    if (!is.vector(k, mode = "numeric") || !all(is.finite(k)) || 
        any(k < 1)) 
        stop(paste("'k' must be a numeric vector", "with all elements greater than or equal to 1", 
            "and no Missing (NA), Infinite (-Inf, Inf),", "or Undefined (Nan) values."))
    if (!is.vector(m, mode = "numeric") || !all(is.finite(m)) || 
        any(m < 1)) 
        stop(paste("'m' must be a numeric vector", "with all elements greater than or equal to 1", 
            "and no Missing (NA), Infinite (Inf, -Inf),", "or Undefined (Nan) values."))
    if (!is.vector(r, mode = "numeric") || !all(is.finite(r)) || 
        any(r < 1)) 
        stop(paste("'r' must be a numeric vector", "with all elements greater than or equal to 1", 
            "and no Missing (NA), Infinite (-Inf, Inf),", "or Undefined (Nan) values."))
    if (pi.type == "upper") 
        lpl.rank <- 0
    else n.plus.one.minus.upl.rank <- 0
    if (!is.vector(lpl.rank, mode = "numeric") || !all(is.finite(lpl.rank)) || 
        any(lpl.rank != trunc(lpl.rank)) || any(lpl.rank < 0 | 
        lpl.rank >= n)) 
        stop(paste("'lpl.rank' must be a vector of non-negative", 
            "integers less than the corresponding value of 'n'"))
    if (pi.type == "lower" & any(lpl.rank < 1)) 
        stop("When pi.type='lower', all values of 'lpl.rank' must be positive integers")
    if (!is.vector(n.plus.one.minus.upl.rank, mode = "numeric") || 
        !all(is.finite(n.plus.one.minus.upl.rank)) || any(n.plus.one.minus.upl.rank != 
        trunc(n.plus.one.minus.upl.rank)) || any(n.plus.one.minus.upl.rank < 
        0 | n.plus.one.minus.upl.rank >= n)) 
        stop(paste("'n.plus.one.minus.upl.rank' must be a vector of non-negative", 
            "integers less than the corresponding value of 'n'"))
    if (pi.type == "upper" & any(n.plus.one.minus.upl.rank < 
        1)) 
        stop("When pi.type='upper' all values of 'n.plus.one.minus.upl.rank' must be positive integers")
    if (!is.vector(delta.over.sigma, mode = "numeric") || any(is.na(delta.over.sigma))) 
        stop(paste("'delta.over.sigma' must be a numeric vector", 
            "with no Missing (NA) or Undefined (Nan) values."))
    if (!is.vector(r.shifted, mode = "numeric") || !all(is.finite(r.shifted)) || 
        !all(r.shifted == trunc(r.shifted)) || any(r.shifted < 
        1) || any(r.shifted > r)) 
        stop(paste("'r.shifted' must be a numeric vector of positive integers", 
            "with all values less than or equal to", "the corresponding values of 'r'"))
    arg.mat <- cbind.no.warn(n = as.vector(n), n.median = as.vector(n.median), 
        k = as.vector(k), m = as.vector(m), r = as.vector(r), 
        lpl.rank = as.vector(lpl.rank), n.plus.one.minus.upl.rank = as.vector(n.plus.one.minus.upl.rank), 
        delta.over.sigma = as.vector(delta.over.sigma), r.shifted = as.vector(r.shifted))
    nrow.arg.mat <- nrow(arg.mat)
    length.rule <- length(rule)
    if (length.rule > nrow.arg.mat) 
        arg.mat <- arg.mat[rep(1:nrow.arg.mat, length.out = length.rule), 
            ]
    else rule <- rep(rule, length.out = nrow.arg.mat)
    for (i in c("n", "n.median", "k", "m", "r", "lpl.rank", "n.plus.one.minus.upl.rank", 
        "delta.over.sigma", "r.shifted")) assign(i, arg.mat[, 
        i])
    index <- rule == "k.of.m"
    if (any(index)) {
        if (any(k[index] > m[index])) 
            stop(paste("For cases where rule='k.of.m',", "all elements of 'k' must be less than or equal to", 
                "the corresponding elements of 'm'"))
    }
    index <- rule == "Modified.CA"
    m[index] <- 4
    N <- length(n)
    conf.level <- numeric(N)
    power <- numeric(N)
    for (i in 1:N) {
        conf.level[i] <- predIntNparSimultaneousConfLevel(n = n[i], 
            n.median = n.median[i], k = k[i], m = m[i], r = r[i], 
            rule = rule[i], lpl.rank = lpl.rank[i], n.plus.one.minus.upl.rank = n.plus.one.minus.upl.rank[i], 
            pi.type = pi.type, integrate.args.list = integrate.args.list)
    }
    if (pi.type == "upper") 
        pl.rank <- n + 1 - n.plus.one.minus.upl.rank
    else pl.rank <- lpl.rank
    if (method == "approx") {
        K <- numeric(N)
        for (i in 1:N) K[i] <- evNormOrdStatsScalar(r = pl.rank[i], 
            n = n[i])
        if (pi.type == "lower") 
            K <- -K
        for (i in 1:N) {
            power[i] <- predIntNormSimultaneousTestPowerScalar(n = n[i], 
                n.mean = n.median[i], K = K[i], k = k[i], m = m[i], 
                r = r.shifted[i], rule = rule[i], delta.over.sigma = delta.over.sigma[i], 
                pi.type = pi.type, conf.level = conf.level[i], 
                integrate.args.list = integrate.args.list)
        }
    }
    else {
        for (i in 1:N) {
            n.i <- n[i]
            n.median.i <- n.median[i]
            k.i <- k[i]
            m.i <- m[i]
            r.i <- r[i]
            pl.rank.i <- pl.rank[i]
            r.shifted.i <- r.shifted[i]
            delta.over.sigma.i <- delta.over.sigma[i]
            mean.i <- c(rep(delta.over.sigma.i, r.shifted.i), 
                rep(0, r.i - r.shifted.i))
            out.vec <- logical(NMC)
            if (pi.type == "upper") {
                test.fcn <- switch(rule[i], k.of.m = function(z, 
                  PL, k) sum(z <= PL) < k, CA = function(z, PL, 
                  k) (z[1] > PL) & (sum(z > PL) >= 2), Modified.CA = function(z, 
                  PL, k) (z[1] > PL) & (sum(z > PL) >= 3))
            }
            else {
                test.fcn <- switch(rule[i], k.of.m = function(z, 
                  PL, k) sum(z >= PL) < k, CA = function(z, PL, 
                  k) (z[1] < PL) & (sum(z < PL) >= 2), Modified.CA = function(z, 
                  PL, k) (z[1] < PL) & (sum(z < PL) >= 3))
            }
            for (j in 1:NMC) {
                x <- rnorm(n.i)
                PL <- sort(x)[pl.rank.i]
                new.x <- array(rnorm(n.median.i * m.i * r.i, 
                  mean = mean.i), dim = c(r.i, n.median.i, m.i))
                new.x <- apply(new.x, c(1, 3), median)
                out.vec[j] <- any(apply(new.x, 1, test.fcn, PL = PL, 
                  k = k.i))
            }
            power[i] <- mean(out.vec)
        }
        if (ci) {
            SE.power <- sqrt(power * (1 - power)/NMC)
            LCL <- power - qnorm(ci.conf.level) * SE.power
            UCL <- power + qnorm(ci.conf.level) * SE.power
            attr(power, "conf.int") <- rbind(LCL = LCL, UCL = UCL)
        }
    }
    power
}
