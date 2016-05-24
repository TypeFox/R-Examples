propTestPower <-
function (n.or.n1, p.or.p1 = 0.5, n2 = n.or.n1, p0.or.p2 = 0.5, 
    alpha = 0.05, sample.type = "one.sample", alternative = "two.sided", 
    approx = TRUE, correct = sample.type == "two.sample", warn = TRUE, 
    return.exact.list = TRUE) 
{
    sample.type <- match.arg(sample.type, c("one.sample", "two.sample"))
    alternative <- match.arg(alternative, c("two.sided", "less", 
        "greater"))
    if (!is.vector(n.or.n1, mode = "numeric") || !is.vector(p.or.p1, 
        mode = "numeric") || !is.vector(p0.or.p2, mode = "numeric") || 
        !is.vector(alpha, mode = "numeric")) 
        stop("'n.or.n1', 'p.or.p1', 'p0.or.p2', and 'alpha' must be numeric vectors.")
    if (!all(is.finite(n.or.n1)) || !all(is.finite(p.or.p1)) || 
        !all(is.finite(p0.or.p2)) || !all(is.finite(alpha))) 
        stop(paste("Missing (NA), Infinite (Inf, -Inf), and", 
            "Undefined (Nan) values are not allowed in", "'n.or.n1', 'p.or.p1', 'p0.or.p2', or 'alpha'"))
    if (any(n.or.n1 < 2)) 
        stop("All values of 'n.or.n1' must be greater than or equal to 2.")
    if (sample.type == "one.sample" && !approx && !all(n.or.n1 == 
        trunc(n.or.n1))) 
        stop(paste("When sample.type='one.sample' and approx=FALSE,", 
            "all values of 'n.or.n1' must be integers"))
    if (any(p.or.p1 <= 0) || any(p.or.p1 >= 1) || any(p0.or.p2 <= 
        0) || any(p0.or.p2 >= 1)) 
        stop("All values of 'p.or.p1' and 'p0.or.p2' must be between 0 and 1")
    if (any(alpha <= 0) || any(alpha >= 1)) 
        stop("All values of 'alpha' must be between 0 and 1.")
    if (sample.type == "two.sample" && !missing(n2)) {
        if (!is.vector(n2, mode = "numeric")) 
            stop("'n2' must be a numeric vector")
        if (!all(is.finite(n2))) 
            stop(paste("Missing (NA), Infinite (Inf, -Inf), and", 
                "Undefined (Nan) values are not allowed in 'n2'"))
        if (any(n2 < 2)) 
            stop("All values of 'n2' must be greater than or equal to 2")
    }
    if (alternative == "two.sided") 
        test.alpha <- alpha/2
    else test.alpha <- alpha
    if (sample.type == "one.sample") {
        arg.mat <- cbind.no.warn(n = as.vector(n.or.n1), p = as.vector(p.or.p1), 
            p0 = as.vector(p0.or.p2), test.alpha = as.vector(test.alpha))
        n <- arg.mat[, "n"]
        p <- arg.mat[, "p"]
        p0 <- arg.mat[, "p0"]
        test.alpha <- arg.mat[, "test.alpha"]
        q0 <- 1 - p0
        q <- 1 - p
        if (approx) {
            if (warn && any(index <- pmin(n * p, n * q) < 5)) 
                warning(paste("The sample size 'n' is too small,", 
                  "relative to the given value of 'p', for the normal", 
                  "approximation to work well for the following", 
                  "element indices:\n\t", paste((1:length(p))[index], 
                    collapse = " "), "\n\t"))
            delta <- p - p0
            cf <- (correct * 1)/(2 * n)
            cf[abs(delta) < cf] <- 0
            z <- qnorm(1 - test.alpha)
            ret.val <- switch(alternative, less = {
                pnorm(-z * sqrt((p0 * q0)/(p * q)) - (delta + 
                  cf)/sqrt((p * q)/n))
            }, greater = {
                1 - pnorm(z * sqrt((p0 * q0)/(p * q)) - (delta - 
                  cf)/sqrt((p * q)/n))
            }, two.sided = {
                pnorm(-z * sqrt((p0 * q0)/(p * q)) - (delta + 
                  cf)/sqrt((p * q)/n)) + 1 - pnorm(z * sqrt((p0 * 
                  q0)/(p * q)) - (delta - cf)/sqrt((p * q)/n))
            })
            names(ret.val) <- NULL
        }
        else {
            switch(alternative, less = {
                q.crit <- qbinom(test.alpha, n, p0)
                index <- pbinom(q.crit, n, p0) > test.alpha
                q.crit[index] <- pmax(q.crit[index] - 1, 0)
                power <- pbinom(q.crit, n, p)
                true.alpha <- pbinom(q.crit, n, p0)
                if (return.exact.list) {
                  names(power) <- names(true.alpha) <- names(q.crit) <- NULL
                  ret.val <- list(power = power, alpha = true.alpha, 
                    q.critical.lower = q.crit)
                } else {
                  names(power) <- NULL
                  ret.val <- power
                }
            }, greater = {
                q.crit <- qbinom(1 - test.alpha, n, p0)
                index <- 1 - pbinom(q.crit, n, p0) > test.alpha
                q.crit[index] <- pmin(q.crit[index] + 1, n[index])
                power <- 1 - pbinom(q.crit, n, p)
                true.alpha <- 1 - pbinom(q.crit, n, p0)
                if (return.exact.list) {
                  names(power) <- names(true.alpha) <- names(q.crit) <- NULL
                  ret.val <- list(power = power, alpha = true.alpha, 
                    q.critical.upper = q.crit)
                } else {
                  names(power) <- NULL
                  ret.val <- power
                }
            }, two.sided = {
                q.crit.lower <- qbinom(test.alpha, n, p0)
                index <- pbinom(q.crit.lower, n, p0) > test.alpha
                q.crit.lower[index] <- pmax(q.crit.lower[index] - 
                  1, 0)
                q.crit.upper <- qbinom(1 - test.alpha, n, p0)
                index <- 1 - pbinom(q.crit.upper, n, p0) > test.alpha
                q.crit.upper[index] <- pmin(q.crit.upper[index] + 
                  1, n[index])
                power <- pbinom(q.crit.lower, n, p) + 1 - pbinom(q.crit.upper, 
                  n, p)
                true.alpha <- pbinom(q.crit.lower, n, p0) + 1 - 
                  pbinom(q.crit.upper, n, p0)
                if (return.exact.list) {
                  names(power) <- names(true.alpha) <- names(q.crit.lower) <- names(q.crit.upper) <- NULL
                  ret.val <- list(power = power, alpha = true.alpha, 
                    q.critical.lower = q.crit.lower, q.critical.upper = q.crit.upper)
                } else {
                  names(power) <- NULL
                  ret.val <- power
                }
            })
            if (any(index <- true.alpha > alpha)) {
                if (return.exact.list) 
                  ret.val$power[index] <- NA
                else ret.val[index] <- NA
                if (warn) 
                  warning(paste("The sample size 'n' is too small", 
                    "relative to the given values of 'alpha' and 'p0'", 
                    "to maintain the significance level for the following", 
                    "element indices:\n\t", paste((1:length(p))[index], 
                      collapse = " "), "\n\t"))
            }
        }
    }
    else {
        arg.mat <- cbind.no.warn(n1 = as.vector(n.or.n1), p1 = as.vector(p.or.p1), 
            n2 = as.vector(n2), p2 = as.vector(p0.or.p2), test.alpha = as.vector(test.alpha))
        n1 <- arg.mat[, "n1"]
        p1 <- arg.mat[, "p1"]
        n2 <- arg.mat[, "n2"]
        p2 <- arg.mat[, "p2"]
        test.alpha <- arg.mat[, "test.alpha"]
        q1 <- 1 - p1
        q2 <- 1 - p2
        if (approx) {
            if (warn && any(index <- pmin(n1 * p1, n1 * q1, n2 * 
                p2, n2 * q2) < 5)) 
                warning(paste("The sample sizes 'n1' and 'n2' are", 
                  "too small, relative to the given values of", 
                  "'p1' and 'p2', for the normal", "approximation to work well for the following", 
                  "element indices:\n\t", paste((1:length(p1))[index], 
                    collapse = " "), "\n\t"))
            delta <- p1 - p2
            n.fac <- 1/n1 + 1/n2
            cf <- correct * 0.5 * n.fac
            cf[abs(delta) < cf] <- 0
            z <- qnorm(1 - test.alpha)
            p.bar <- (n1 * p1 + n2 * p2)/(n1 + n2)
            q.bar <- 1 - p.bar
            ret.val <- switch(alternative, less = {
                pnorm((-z * sqrt(p.bar * q.bar * n.fac) - delta - 
                  cf)/sqrt((p1 * q1)/n1 + (p2 * q2)/n2))
            }, greater = {
                1 - pnorm((z * sqrt(p.bar * q.bar * n.fac) - 
                  delta + cf)/sqrt((p1 * q1)/n1 + (p2 * q2)/n2))
            }, two.sided = {
                pnorm((-z * sqrt(p.bar * q.bar * n.fac) - delta - 
                  cf)/sqrt((p1 * q1)/n1 + (p2 * q2)/n2)) + 1 - 
                  pnorm((z * sqrt(p.bar * q.bar * n.fac) - delta + 
                    cf)/sqrt((p1 * q1)/n1 + (p2 * q2)/n2))
            })
            names(ret.val) <- NULL
        }
        else {
            stop("Power of exact method not available for the two-sample test")
        }
    }
    ret.val
}
