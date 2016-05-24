ciBinomHalfWidth <-
function (n.or.n1, p.hat.or.p1.hat = 0.5, n2 = n.or.n1, p2.hat = 0.4, 
    conf.level = 0.95, sample.type = "one.sample", ci.method = "score", 
    correct = TRUE, warn = TRUE) 
{
    if (missing(sample.type)) {
        sample.type <- ifelse(!missing(n2) || !missing(p2.hat), 
            "two.sample", "one.sample")
    }
    sample.type <- match.arg(sample.type, c("one.sample", "two.sample"))
    ci.method <- match.arg(ci.method, c("score", "exact", "adjusted Wald", 
        "Wald"))
    if (!is.vector(n.or.n1, mode = "numeric") || !is.vector(p.hat.or.p1.hat, 
        mode = "numeric") || !is.vector(conf.level, mode = "numeric")) 
        stop("'n.or.n1', 'p.hat', and 'conf.level' must be numeric vectors.")
    if (!all(is.finite(n.or.n1)) || !all(is.finite(p.hat.or.p1.hat)) || 
        !all(is.finite(conf.level))) 
        stop(paste("Missing (NA), Infinite (Inf, -Inf), and", 
            "Undefined (Nan) values are not allowed in", "'n.or.n1', 'p.hat.or.p1.hat', or 'conf.level'"))
    if (any(n.or.n1 < 2)) 
        stop("All values of 'n.or.n1' must be greater than or equal to 2")
    if (any(p.hat.or.p1.hat < .Machine$double.eps) || any(p.hat.or.p1.hat > 
        1 - .Machine$double.eps)) 
        stop("All values of 'p.hat.or.p1.hat' must be strictly between 0 and 1.")
    if (any(conf.level <= .Machine$double.eps) || any(conf.level >= 
        1 - .Machine$double.eps)) 
        stop("All values of 'conf.level' must be greater than 0 and less than 1")
    if (sample.type == "two.sample") {
        if (ci.method == "exact") 
            stop("Exact method not available for two-sample confidence interval")
        if (!missing(n2)) {
            if (!is.vector(n2, mode = "numeric")) 
                stop("'n2' must be a numeric vector")
            if (!all(is.finite(n2))) 
                stop(paste("Missing (NA), Infinite (Inf, -Inf), and", 
                  "Undefined (Nan) values are not allowed in 'n2'"))
            if (any(n2 < 2)) 
                stop("All values of 'n2' must be greater than or equal to 2")
        }
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
        arg.mat <- cbind.no.warn(n1 = as.vector(n.or.n1), n2 = as.vector(n2), 
            p1.hat = as.vector(p.hat.or.p1.hat), p2.hat = as.vector(p2.hat), 
            conf.level = as.vector(conf.level))
        for (i in c("n1", "n2", "p1.hat", "p2.hat", "conf.level")) assign(i, 
            arg.mat[, i])
        p1.hat.high <- ceiling(n1 * p1.hat)/n1
        p1.hat.low <- floor(n1 * p1.hat)/n1
        index <- p1.hat.high != p1.hat.low
        index.tied <- abs((p1.hat.high - p1.hat) - (p1.hat - 
            p1.hat.low)) <= .Machine$double.eps
        index.high <- (index & !index.tied & ((p1.hat.high - 
            p1.hat) < (p1.hat - p1.hat.low))) | (index & index.tied & 
            (abs(p1.hat.high - 0.5) <= abs(p1.hat.low - 0.5)))
        index.low <- (index & !index.tied & ((p1.hat.high - p1.hat) > 
            (p1.hat - p1.hat.low))) | (index & index.tied & (abs(p1.hat.high - 
            0.5) > abs(p1.hat.low - 0.5)))
        p1.hat[index.high] <- p1.hat.high[index.high]
        p1.hat[index.low] <- p1.hat.low[index.low]
        x1 <- round(n1 * p1.hat)
        p2.hat.high <- ceiling(n2 * p2.hat)/n2
        p2.hat.low <- floor(n2 * p2.hat)/n2
        index <- p2.hat.high != p2.hat.low
        index.tied <- abs((p2.hat.high - p2.hat) - (p2.hat - 
            p2.hat.low)) <= .Machine$double.eps
        index.high <- (index & !index.tied & ((p2.hat.high - 
            p2.hat) < (p2.hat - p2.hat.low))) | (index & index.tied & 
            (abs(p2.hat.high - 0.5) <= abs(p2.hat.low - 0.5)))
        index.low <- (index & !index.tied & ((p2.hat.high - p2.hat) > 
            (p2.hat - p2.hat.low))) | (index & index.tied & (abs(p2.hat.high - 
            0.5) > abs(p2.hat.low - 0.5)))
        p2.hat[index.high] <- p2.hat.high[index.high]
        p2.hat[index.low] <- p2.hat.low[index.low]
        x2 <- round(n2 * p2.hat)
        switch(ci.method, score = {
            N <- length(n1)
            half.width <- numeric(N)
            o.warn <- options(warn = -1)
            for (i in 1:N) {
                x1.i <- x1[i]
                n1.i <- n1[i]
                x2.i <- x2[i]
                n2.i <- n2[i]
                conf.i <- conf.level[i]
                half.width[i] <- diff(prop.test(x = c(x1.i, x2.i), 
                  n = c(n1.i, n2.i), conf.level = conf.i, correct = correct)$conf.int)/2
            }
            options(o.warn)
        }, `adjusted Wald` = {
            const <- qnorm(1 - (1 - conf.level)/2)^2
            x1.a <- x1 + const/4
            x2.a <- x2 + const/4
            n1.a <- n1 + const/2
            n2.a <- n2 + const/2
            p1.hat.a <- x1.a/n1.a
            p2.hat.a <- x2.a/n2.a
            q1.hat.a <- 1 - p1.hat.a
            q2.hat.a <- 1 - p2.hat.a
            s.hat <- sqrt((p1.hat.a * q1.hat.a)/n1.a + (p2.hat.a * 
                q2.hat.a)/n2.a)
            half.width <- qnorm(1 - (1 - conf.level)/2) * s.hat
        }, Wald = {
            q1.hat <- 1 - p1.hat
            q2.hat <- 1 - p2.hat
            s.hat <- sqrt((p1.hat * q1.hat)/n1 + (p2.hat * q2.hat)/n2)
            cf <- correct * 0.5 * (1/n1 + 1/n2)
            half.width <- qnorm(1 - (1 - conf.level)/2) * s.hat + 
                cf
            if (warn && any(index <- pmin(n1 * p1.hat, n1 * q1.hat, 
                n2 * p2.hat, n2 * q2.hat) < 5)) warning(paste("The sample sizes 'n1' and 'n2' are too small,", 
                "relative to the given values of 'p1.hat' and 'p2.hat',", 
                "for the normal approximation to work well for the following", 
                "element indices:\n\t", paste((1:length(half.width))[index], 
                  collapse = " "), "\n\t"))
        })
        names(half.width) <- names(n1) <- names(n2) <- names(p1.hat) <- names(p2.hat) <- NULL
        ret.val <- list(half.width = half.width, n1 = n1, p1.hat = p1.hat, 
            n2 = n2, p2.hat = p2.hat)
    }
    else {
        arg.mat <- cbind.no.warn(n = as.vector(n.or.n1), p.hat = as.vector(p.hat.or.p1.hat), 
            conf.level = as.vector(conf.level))
        for (i in c("n", "p.hat", "conf.level")) assign(i, arg.mat[, 
            i])
        p.hat.high <- ceiling(n * p.hat)/n
        p.hat.low <- floor(n * p.hat)/n
        index <- p.hat.high != p.hat.low
        index.tied <- abs((p.hat.high - p.hat) - (p.hat - p.hat.low)) <= 
            .Machine$double.eps
        index.high <- (index & !index.tied & ((p.hat.high - p.hat) < 
            (p.hat - p.hat.low))) | (index & index.tied & (abs(p.hat.high - 
            0.5) <= abs(p.hat.low - 0.5)))
        index.low <- (index & !index.tied & ((p.hat.high - p.hat) > 
            (p.hat - p.hat.low))) | (index & index.tied & (abs(p.hat.high - 
            0.5) > abs(p.hat.low - 0.5)))
        p.hat[index.high] <- p.hat.high[index.high]
        p.hat[index.low] <- p.hat.low[index.low]
        x <- round(n * p.hat)
        N <- length(n)
        half.width <- numeric(N)
        switch(ci.method, score = {
            for (i in 1:N) {
                x.i <- x[i]
                n.i <- n[i]
                conf.i <- conf.level[i]
                half.width[i] <- diff(ebinom(x = x.i, size = n.i, 
                  ci = TRUE, ci.type = "two-sided", ci.method = "score", 
                  conf.level = conf.i, correct = correct)$interval$limits)/2
            }
        }, exact = {
            for (i in 1:N) {
                x.i <- x[i]
                n.i <- n[i]
                conf.i <- conf.level[i]
                half.width[i] <- diff(ebinom(x = x.i, size = n.i, 
                  ci = TRUE, ci.type = "two-sided", ci.method = "exact", 
                  conf.level = conf.i)$interval$limits)/2
            }
        }, `adjusted Wald` = {
            for (i in 1:N) {
                x.i <- x[i]
                n.i <- n[i]
                conf.i <- conf.level[i]
                half.width[i] <- diff(ebinom(x = x.i, size = n.i, 
                  ci = TRUE, ci.type = "two-sided", ci.method = "adjusted Wald", 
                  conf.level = conf.i)$interval$limits)/2
            }
        }, Wald = {
            cf <- (correct * 0.5)/n
            half.width <- qnorm(1 - (1 - conf.level)/2) * sqrt((p.hat * 
                (1 - p.hat))/n) + cf
            if (warn && any(index <- pmin(n * p.hat, n * (1 - 
                p.hat)) < 5)) warning(paste("The sample size 'n' is too small,", 
                "relative to the given value of 'p.hat', for the normal", 
                "approximation to work well for the following", 
                "element indices:\n\t", paste((1:length(half.width))[index], 
                  collapse = " "), "\n\t"))
        })
        names(half.width) <- names(n) <- names(p.hat) <- NULL
        ret.val <- list(half.width = half.width, n = n, p.hat = p.hat)
    }
    method <- switch(ci.method, score = {
        ifelse(!correct, "Score normal approximation", "Score normal approximation, with continuity correction")
    }, exact = "Exact", Wald = {
        ifelse(!correct, "Wald normal approximation", "Wald normal approximation, with continuity correction")
    }, `adjusted Wald` = "Adjusted Wald normal approximation")
    c(ret.val, method = method)
}
