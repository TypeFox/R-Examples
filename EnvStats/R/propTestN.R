propTestN <-
function (p.or.p1, p0.or.p2, alpha = 0.05, power = 0.95, sample.type = "one.sample", 
    alternative = "two.sided", ratio = 1, approx = TRUE, correct = sample.type == 
        "two.sample", round.up = TRUE, warn = TRUE, return.exact.list = TRUE, 
    n.min = 2, n.max = 10000, tol.alpha = 0.1 * alpha, tol = 1e-07, 
    maxiter = 1000) 
{
    sample.type <- match.arg(sample.type, c("one.sample", "two.sample"))
    alternative <- match.arg(alternative, c("two.sided", "less", 
        "greater"))
    if (!is.vector(p.or.p1, mode = "numeric") || !is.vector(alpha, 
        mode = "numeric") || !is.vector(power, mode = "numeric") || 
        !is.vector(ratio, mode = "numeric")) 
        stop("'p.or.p1', 'alpha', 'power', and 'ratio' must be numeric vectors")
    if (!all(is.finite(p.or.p1)) || !all(is.finite(alpha)) || 
        !all(is.finite(power)) || !all(is.finite(ratio))) 
        stop(paste("Missing (NA), Infinite (Inf, -Inf), and", 
            "Undefined (Nan) values are not allowed in", "'p.or.p1', 'alpha', 'power', or 'ratio'"))
    if (any(p.or.p1 <= 0) || any(p.or.p1 >= 1)) 
        stop("All values of 'p.or.p1' must be greater than 0 and less than 1")
    if (any(alpha <= 0) || any(alpha >= 1)) 
        stop("All values of 'alpha' must be greater than 0 and less than 1")
    if (any(power <= alpha) || any(power >= 1)) 
        stop(paste("All values of 'power' must be greater than", 
            "the corresponding elements of 'alpha' and less than 1"))
    if (alternative == "greater" && any(p.or.p1 <= p0.or.p2)) 
        stop(paste("When alternative='greater', all values of 'p.or.p1'", 
            "must be greater than the corresponding values of 'p0.or.p2'"))
    if (alternative == "less" && any(p.or.p1 >= p0.or.p2)) 
        stop(paste("When alternative='less', all values of 'p.or.p1'", 
            "must be less than the corresponding values of 'p0.or.p2'"))
    if (alternative == "two.sided" && any(p.or.p1 == p0.or.p2)) 
        stop(paste("When alternative='two.sided', all values of 'p.or.p1'", 
            "must be unequal to the corresponding values of 'p0.or.p2'"))
    if (length(maxiter) != 1 || !is.numeric(maxiter) || !is.finite(maxiter) || 
        maxiter < 2 || maxiter != round(maxiter)) 
        stop("'maxiter' must be an integer greater than 1")
    if (sample.type == "two.sample") {
        if (!approx) 
            stop("Sample size based on exact power for two-sample test not available")
        if (!is.vector(ratio, mode = "numeric")) 
            stop("'ratio' must be a numeric vector.")
        if (!all(is.finite(ratio))) 
            stop(paste("Missing (NA), Infinite (Inf, -Inf), and", 
                "Undefined (Nan) values are not allowed in", 
                "'ratio'"))
        if (ratio.gt.1 <- !all(ratio == 1)) {
            if (any(ratio < 1)) 
                stop("All values of 'ratio' must be greater than or equal to 1")
        }
        arg.mat <- cbind.no.warn(power = as.vector(power), p1 = as.vector(p.or.p1), 
            p2 = as.vector(p0.or.p2), alpha = as.vector(alpha), 
            r = as.vector(ratio))
        power <- arg.mat[, "power"]
        p1 <- arg.mat[, "p1"]
        p2 <- arg.mat[, "p2"]
        alpha <- arg.mat[, "alpha"]
        r <- arg.mat[, "r"]
        if (alternative == "two.sided") 
            test.alpha <- alpha/2
        else test.alpha <- alpha
        p.bar <- (p1 + r * p2)/(r + 1)
        q.bar <- 1 - p.bar
        delta <- abs(p1 - p2)
        n.vec <- (qnorm(1 - test.alpha) * sqrt((r + 1) * p.bar * 
            q.bar) + qnorm(power) * sqrt(r * p1 * (1 - p1) + 
            p2 * (1 - p2)))^2/(r * delta^2)
        if (correct) 
            n.vec <- (n.vec/4) * (1 + sqrt(1 + (2 * (r + 1))/(r * 
                n.vec * delta)))^2
        n2 <- ratio * n.vec
        if (round.up) {
            n.vec <- ceiling(n.vec)
            n2 <- ceiling(n2)
        }
        if (warn && any(index <- pmin(n.vec * p1, n.vec * (1 - 
            p1), n2 * p2, n2 * (1 - p2)) < 5)) 
            warning(paste("The computed sample sizes 'n1' and 'n2'", 
                "are too small, relative to the given values of", 
                "'p1' and 'p2', for the normal", "approximation to work well for the following", 
                "element indices:\n\t", paste((1:length(p1))[index], 
                  collapse = " "), "\n\t"))
        if (ratio.gt.1) {
            names(n.vec) <- names(n2) <- NULL
            ret.val <- list(n1 = n.vec, n2 = n2)
        }
        else {
            names(n.vec) <- NULL
            ret.val <- n.vec
        }
    }
    else {
        arg.mat <- cbind.no.warn(power = as.vector(power), p0 = as.vector(p0.or.p2), 
            p = as.vector(p.or.p1), alpha = as.vector(alpha), 
            tol.alpha = as.vector(tol.alpha))
        power <- arg.mat[, "power"]
        p0 <- arg.mat[, "p0"]
        p <- arg.mat[, "p"]
        alpha <- arg.mat[, "alpha"]
        tol.alpha <- arg.mat[, "tol.alpha"]
        if (approx) {
            if (alternative == "two.sided") 
                test.alpha <- alpha/2
            else test.alpha <- alpha
            q0 <- 1 - p0
            q <- 1 - p
            n.vec <- (qnorm(power) * sqrt(p * q) - qnorm(test.alpha) * 
                sqrt(p0 * q0))^2/((p0 - p)^2)
            if (correct) 
                n.vec <- (n.vec/4) * (1 + sqrt(1 + 2/(n.vec * 
                  abs(p0 - p))))^2
            ceiling.n.vec <- ceiling(n.vec)
            if (round.up) 
                n.vec <- ceiling.n.vec
            names(n.vec) <- NULL
            if (warn && any(index <- pmin(n.vec * p, n.vec * 
                q) < 5)) 
                warning(paste("The computed sample size 'n' is too small,", 
                  "relative to the given value of 'p', for the normal", 
                  "approximation to work well for the following", 
                  "element indices:\n\t", paste((1:length(p))[index], 
                    collapse = " "), "\n\t"))
            ret.val <- n.vec
        }
        else {
            if (length(n.min) != 1 || length(n.max) != 1 || !is.numeric(n.min) || 
                !is.numeric(n.max) || !is.finite(n.min) || !is.finite(n.max) || 
                n.min < 2 || n.max <= n.min || n.min != round(n.min) || 
                n.max != round(n.max)) 
                stop(paste("'n.min' must be an integer greater than 1", 
                  "and 'n.max' must be an integer greater than 'n.min'"))
            if (any(tol.alpha <= 0) || any(tol.alpha >= alpha)) 
                stop(paste("All values of 'tol.alpha' must be greater than 0 and", 
                  "less than the corresponding value of 'alpha'."))
            if (any(alpha >= 0.5) && any(tol.alpha >= (1 - alpha))) 
                stop(paste("When an element of 'alpha' is greater than", 
                  "or equal to 0.5, the corresponding element of 'tol.alpha'", 
                  "must be less than 1 - 'alpha'"))
            N <- length(power)
            dum.list <- propTestPower(n.or.n1 = n.max, p.or.p1 = p, 
                p0.or.p2 = p0, alpha = alpha, sample.type = "one.sample", 
                alternative = alternative, approx = FALSE, return.exact.list = TRUE)
            dum.power <- dum.list$power
            dum.alpha <- dum.list$alpha
            index <- (dum.power < power) | (abs(dum.alpha - alpha) > 
                tol.alpha)
            if (any(index)) 
                stop(paste("The supplied value of 'n.max' is too small,", 
                  "relative to the given values of 'power', 'alpha', and 'tol.alpha',", 
                  "in order to achieve the desired power for the following", 
                  "element indices:\n\t", paste((1:N)[index], 
                    collapse = " "), "\n\t"))
            results.mat <- matrix(as.numeric(NA), nrow = N, ncol = 5)
            dimnames(results.mat) <- list(NULL, c("n", "power", 
                "alpha", "q.critical.lower", "q.critical.upper"))
            if (alternative == "greater") 
                results.mat <- results.mat[, c("n", "power", 
                  "alpha", "q.critical.upper"), drop = FALSE]
            if (alternative == "less") 
                results.mat <- results.mat[, c("n", "power", 
                  "alpha", "q.critical.lower"), drop = FALSE]
            for (i in 1:N) {
                results.mat[i, ] <- propTestNVecExact(p = p[i], 
                  p0 = p0[i], alpha = alpha[i], power = power[i], 
                  alternative = alternative, round.up = round.up, 
                  tol = tol, n.min = n.min, n.max = n.max, tol.alpha = tol.alpha[i], 
                  maxiter = maxiter)
            }
            if (return.exact.list) {
                ret.val <- unclass(data.frame(results.mat))
                attr(ret.val, "row.names") <- NULL
            }
            else ret.val <- results.mat[, "n"]
        }
    }
    ret.val
}
