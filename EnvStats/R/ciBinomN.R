ciBinomN <-
function (half.width, p.hat.or.p1.hat = 0.5, p2.hat = 0.4, conf.level = 0.95, 
    sample.type = "one.sample", ratio = 1, ci.method = "score", 
    correct = TRUE, warn = TRUE, n.or.n1.min = 2, n.or.n1.max = 10000, 
    tol.half.width = 5e-04, tol.p.hat = 5e-04, tol = 1e-07, maxiter = 1000) 
{
    if (missing(sample.type)) {
        sample.type <- ifelse(!missing(p2.hat) || !missing(ratio), 
            "two.sample", "one.sample")
    }
    sample.type <- match.arg(sample.type, c("one.sample", "two.sample"))
    ci.method <- match.arg(ci.method, c("score", "exact", "adjusted Wald", 
        "Wald"))
    if (!is.vector(half.width, mode = "numeric") || !is.vector(p.hat.or.p1.hat, 
        mode = "numeric") || !is.vector(conf.level, mode = "numeric")) 
        stop("'half.width', 'p.hat.or.p1.hat', and 'conf.level' must be numeric vectors.")
    if (!all(is.finite(half.width)) || !all(is.finite(p.hat.or.p1.hat)) || 
        !all(is.finite(conf.level))) 
        stop(paste("Missing (NA), Infinite (Inf, -Inf), and", 
            "Undefined (Nan) values are not allowed in", "'half.width', 'p.hat.or.p1.hat', or 'conf.level'"))
    if (any(half.width < .Machine$double.eps) || any(half.width >= 
        0.5)) 
        stop(paste("All values of 'half.width' must be positive", 
            "and less than 0.5"))
    if (any(p.hat.or.p1.hat < .Machine$double.eps) || any(p.hat.or.p1.hat > 
        1 - .Machine$double.eps)) 
        stop("All values of 'p.hat.or.p1.hat' must be strictly between 0 and 1.")
    if (any(conf.level <= .Machine$double.eps) || any(conf.level >= 
        1 - .Machine$double.eps)) 
        stop("All values of 'conf.level' must be greater than 0 and less than 1")
    if (length(n.or.n1.min) != 1 || length(n.or.n1.max) != 1 || 
        !is.numeric(n.or.n1.min) || !is.numeric(n.or.n1.max) || 
        !is.finite(n.or.n1.min) || !is.finite(n.or.n1.max) || 
        n.or.n1.min < 2 || n.or.n1.max <= n.or.n1.min || n.or.n1.min != 
        round(n.or.n1.min) || n.or.n1.max != round(n.or.n1.max)) 
        stop(paste("'n.or.n1.min' must be an integer greater than 1", 
            "and 'n.or.n1.max' must be an integer greater than 'n.or.n1.min'"))
    if (length(tol.half.width) != 1 || tol.half.width <= .Machine$double.eps || 
        tol.half.width >= 0.5) 
        stop("'tol.half.width' must be a scalar greater than 0 and less than 0.5")
    if (length(tol.p.hat) != 1 || tol.p.hat <= .Machine$double.eps || 
        tol.p.hat >= 0.5) 
        stop("'tol.p.hat' must be a scalar greater than 0 and less than 0.5")
    if (sample.type == "two.sample") {
        if (ci.method == "exact") 
            stop("Exact method not available for two-sample confidence interval")
        if (!missing(p2.hat)) {
            if (!is.vector(p2.hat, mode = "numeric")) 
                stop("'p2.hat' must be a numeric vector")
            if (!all(is.finite(p2.hat))) 
                stop(paste("Missing (NA), Infinite (Inf, -Inf), and", 
                  "Undefined (Nan) values are not allowed in 'p2.hat'"))
            if (any(p2.hat < .Machine$double.eps) || any(p2.hat > 
                1 - .Machine$double.eps)) 
                stop("All values of 'p2.hat' must be strictly between 0 and 1")
        }
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
        arg.mat <- cbind.no.warn(half.width = as.vector(half.width), 
            p1.hat = as.vector(p.hat.or.p1.hat), p2.hat = as.vector(p2.hat), 
            conf.level = as.vector(conf.level), ratio = as.vector(ratio))
        half.width <- arg.mat[, "half.width"]
        p1.hat <- arg.mat[, "p1.hat"]
        p2.hat <- arg.mat[, "p2.hat"]
        conf.level <- arg.mat[, "conf.level"]
        ratio <- arg.mat[, "ratio"]
        if (any(((p1.hat + tol.p.hat) - (p2.hat - tol.p.hat) + 
            half.width + tol.half.width) > 1) || any(((p1.hat - 
            tol.p.hat) - (p2.hat + tol.p.hat) - half.width - 
            tol.half.width) < -1)) {
            stop(paste("The values of \"p.hat.or.p1.hat\", \"p2.hat\", \"tol.p.hat\",", 
                "\"half.width\", and \"tol.half.width\" must satisfy\n\t", 
                "((p.hat.or.p1.hat + tol.p.hat) - (p2.hat - tol.p.hat) + half.width + tol.half.width) <= 1, and\n\t", 
                "((p.hat.or.p1.hat - tol.p.hat) - (p2.hat + tol.p.hat) - half.width - tol.half.width) >= -1."))
        }
        N <- nrow(arg.mat)
        dum.list <- ciBinomHalfWidth(n.or.n1 = n.or.n1.max, p.hat.or.p1.hat = p1.hat, 
            n2 = round(ratio * n.or.n1.max), p2.hat = p2.hat, 
            conf.level = conf.level, sample.type = "two.sample", 
            ci.method = ci.method, correct = correct, warn = FALSE)
        index <- dum.list$half.width <= (half.width + tol.half.width) & 
            abs(p1.hat - dum.list$p1.hat) <= tol.p.hat & abs(p2.hat - 
            dum.list$p2.hat) <= tol.p.hat
        if (!all(index)) {
            if (N == 1) 
                stop(paste("The maximum allowed sample size 'n.or.n1.max' is\n", 
                  "too small to achieve the desired half-width for the given values\n", 
                  "of 'tol.half.width' and 'tol.p.hat'.  Try increasing the values of\n", 
                  "'n.or.n1.max', 'half.width', 'tol.half.width' and/or 'tol.p.hat'."))
            else stop(paste("The maximum allowed sample size 'n.or.n1.max' is\n", 
                "too small to achieve the desired half-width for the given values\n", 
                "of 'tol.half.width' and 'tol.p.hat' for the following indices:\n", 
                paste((1:N)[!index], collapse = " "), ".\nTry increasing the values of\n", 
                "'n.or.n1.max', 'half.width', 'tol.half.width' and/or 'tol.p.hat'.", 
                sep = ""))
        }
        result.mat <- matrix(as.numeric(NA), nrow = N, ncol = 5)
        dimnames(result.mat) <- list(NULL, c("n1", "p1.hat", 
            "n2", "p2.hat", "half.width"))
        for (i in 1:N) {
            result.mat[i, ] <- ciBinomN.vec(half.width = half.width[i], 
                p.hat.or.p1.hat = p1.hat[i], p2.hat = p2.hat[i], 
                conf.level = conf.level[i], sample.type = "two.sample", 
                ratio = ratio[i], ci.method = ci.method, correct = correct, 
                n.or.n1.min = n.or.n1.min, n.or.n1.max = n.or.n1.max, 
                tol.half.width = tol.half.width, tol.p.hat = tol.p.hat, 
                tol = tol, maxiter = maxiter)
        }
        n1 <- result.mat[, "n1"]
        p1.hat <- result.mat[, "p1.hat"]
        n2 <- result.mat[, "n2"]
        p2.hat <- result.mat[, "p2.hat"]
        half.width <- result.mat[, "half.width"]
        if (ci.method == "Wald" && warn && any(index <- pmin(n1 * 
            p1.hat, n1 * (1 - p1.hat), n2 * p2.hat, n2 * (1 - 
            p2.hat)) < 5)) 
            warning(paste("The computed sample size is too small, relative to", 
                "the values of 'p1.hat' and 'p2.hat', for the normal approximation", 
                "to work well for the following element indices:\n\t", 
                paste((1:N)[index], collapse = " "), "\n\t"))
        ret.val <- list(n1 = as.vector(n1), n2 = as.vector(n2), 
            p1.hat = as.vector(p1.hat), p2.hat = as.vector(p2.hat), 
            half.width = as.vector(half.width))
    }
    else {
        arg.mat <- cbind.no.warn(half.width = as.vector(half.width), 
            p.hat = as.vector(p.hat.or.p1.hat), conf.level = as.vector(conf.level))
        half.width <- arg.mat[, "half.width"]
        p.hat <- arg.mat[, "p.hat"]
        conf.level <- arg.mat[, "conf.level"]
        if (any((p.hat + tol.p.hat + half.width + tol.half.width) > 
            1) || any((p.hat - tol.p.hat - half.width - tol.half.width) < 
            0)) {
            stop(paste("The values of \"p.hat.or.p1.hat\", \"tol.p.hat\",", 
                "\"half.width\", and \"tol.half.width\" must satisfy\n\t", 
                "(p.hat.or.p1.hat + tol.p.hat + half.width + tol.half.width) <= 1, and\n\t", 
                "(p.hat.or.p1.hat - tol.p.hat - half.width - tol.half.width) >= 0."))
        }
        N <- nrow(arg.mat)
        dum.list <- ciBinomHalfWidth(n.or.n1 = n.or.n1.max, p.hat.or.p1.hat = p.hat, 
            conf.level = conf.level, sample.type = "one.sample", 
            ci.method = ci.method, correct = correct, warn = FALSE)
        index <- dum.list$half.width <= (half.width + tol.half.width) & 
            abs(p.hat - dum.list$p.hat) <= tol.p.hat
        if (!all(index)) {
            if (N == 1) 
                stop(paste("The maximum allowed sample size 'n.or.n1.max' is\n", 
                  "too small to achieve the desired half-width for the given values\n", 
                  "of 'tol.half.width' and 'tol.p.hat'.  Try increasing the values of\n", 
                  "'n.or.n1.max', 'half.width', 'tol.half.width' and/or 'tol.p.hat'."))
            else stop(paste("The maximum allowed sample size 'n.or.n1.max' is\n", 
                "too small to achieve the desired half-width for the given values\n", 
                "of 'tol.half.width' and 'tol.p.hat' for the following indices:\n", 
                paste((1:N)[!index], collapse = " "), ".\nTry increasing the values of\n", 
                "'n.or.n1.max', 'half.width', 'tol.half.width' and/or 'tol.p.hat'.", 
                sep = "", collapse = ""))
        }
        result.mat <- matrix(as.numeric(NA), nrow = N, ncol = 3)
        dimnames(result.mat) <- list(NULL, c("n", "p.hat", "half.width"))
        for (i in 1:N) {
            result.mat[i, ] <- ciBinomN.vec(half.width = half.width[i], 
                p.hat.or.p1.hat = p.hat[i], conf.level = conf.level[i], 
                sample.type = "one.sample", ci.method = ci.method, 
                correct = correct, n.or.n1.min = n.or.n1.min, 
                n.or.n1.max = n.or.n1.max, tol.half.width = tol.half.width, 
                tol.p.hat = tol.p.hat, tol = tol, maxiter = maxiter)
        }
        n <- result.mat[, "n"]
        p.hat <- result.mat[, "p.hat"]
        half.width <- result.mat[, "half.width"]
        if (ci.method == "Wald" && warn && any(index <- pmin(n * 
            p.hat, n * (1 - p.hat)) < 5)) 
            warning(paste("The computed sample size is too small, relative to", 
                "the value of 'p.hat', for the normal approximation", 
                "to work well for the following element indices:\n\t", 
                paste((1:N)[index], collapse = " "), "\n\t"))
        ret.val <- list(n = as.vector(n), p.hat = as.vector(p.hat), 
            half.width = as.vector(half.width))
    }
    method <- switch(ci.method, score = {
        ifelse(!correct, "Score normal approximation", "Score normal approximation, with continuity correction")
    }, exact = "Exact", Wald = {
        ifelse(!correct, "Wald normal approximation", "Wald normal approximation, with continuity correction")
    }, `adjusted Wald` = "Adjusted Wald normal approximation")
    c(ret.val, list(method = method))
}
