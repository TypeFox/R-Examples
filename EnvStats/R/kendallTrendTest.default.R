kendallTrendTest.default <-
function (y, x = seq(along = y), alternative = "two.sided", correct = TRUE, 
    ci.slope = TRUE, conf.level = 0.95, warn = TRUE, data.name = NULL, 
    data.name.x = NULL, parent.of.data = NULL, subset.expression = NULL, 
    ...) 
{
    if (is.null(data.name)) 
        data.name <- deparse(substitute(y))
    y <- as.vector(unlist(y))
    if (!is.numeric(y)) 
        stop("All elements of 'y' must be a numeric")
    miss <- FALSE
    if (missing(x)) {
        data.name.x <- NULL
        if ((bad.obs <- sum(!(y.ok <- is.finite(y)))) > 0) {
            if (warn) {
                is.not.finite.warning(y)
                warning(paste(bad.obs, "observations with NA/NaN/Inf in 'y' removed."))
            }
            x <- x[y.ok]
            y <- y[y.ok]
        }
        n <- length(y)
        if (n < 2) {
            if (warn) 
                warning(paste("'y' does not contain at least two non-missing,", 
                  "finite observations"))
            miss <- TRUE
        }
    }
    else {
        if (is.null(data.name.x)) 
            data.name.x <- deparse(substitute(x))
        x <- as.vector(unlist(x))
        if (!is.numeric(x) || length(x) != length(y)) 
            stop(paste("All elements of 'x' must be a numeric, and", 
                "'x' must have the same number of elements as 'y'"))
        names(data.name.x) <- "x"
        names(data.name) <- "y"
        if ((bad.obs <- sum(!(ok <- is.finite(x) & is.finite(y)))) > 
            0) {
            if (warn) {
                is.not.finite.warning(x)
                is.not.finite.warning(y)
                warning(paste(bad.obs, "observations with NA/NaN/Inf in 'x' and 'y' removed."))
            }
            x <- x[ok]
            y <- y[ok]
        }
        n <- length(y)
        if (n < 2) {
            if (warn) 
                warning(paste("'x' and 'y' do not contain at least two non-missing,", 
                  "finite observations"))
            miss <- TRUE
        }
        else if (length(unique(x)) < 2) {
            if (warn) 
                warning("'x' does not contain at least 2 distinct values")
            miss <- TRUE
        }
    }
    alternative <- match.arg(alternative, c("two.sided", "greater", 
        "less"))
    if (ci.slope && !is.vector(conf.level, mode = "numeric") || 
        length(conf.level) != 1 || conf.level <= 0 || conf.level >= 
        1) 
        stop("'conf.level' must be a numeric scalar between 0 and 1")
    if (miss) {
        stat <- c(z = NA)
        p.value <- NA
        estimate <- c(tau = NA, slope = NA, intercept = NA)
        method <- "Kendall's Test for Trend"
        estimation.method <- paste("slope:      Theil/Sen Estimator", 
            "intercept:  Conover's Estimator", sep = paste("\n", 
                space(33), sep = ""))
        S <- NA
        var.S <- NA
        slopes <- NA
    }
    else {
        vark <- function(x, y) {
            ties.x <- rle(sort(x))$lengths
            ties.y <- rle(sort(y))$lengths
            n <- length(x)
            t1 <- n * (n - 1) * (2 * n + 5)
            t2 <- sum(ties.x * (ties.x - 1) * (2 * ties.x + 5))
            t3 <- sum(ties.y * (ties.y - 1) * (2 * ties.y + 5))
            v1 <- (t1 - t2 - t3)/18
            if (n > 2) {
                t1 <- sum(ties.x * (ties.x - 1) * (ties.x - 2))
                t2 <- sum(ties.y * (ties.y - 1) * (ties.y - 2))
                v2 <- (t1 * t2)/(9 * n * (n - 1) * (n - 2))
            }
            else v2 <- 0
            t1 <- sum(ties.x * (ties.x - 1)) * sum(ties.y * (ties.y - 
                1))
            v3 <- t1/(2 * n * (n - 1))
            v1 + v2 + v3
        }
        index <- 2:n
        S <- sum(sapply(index, function(i, x, y) {
            sum(sign((x[i] - x[1:(i - 1)]) * (y[i] - y[1:(i - 
                1)])))
        }, x, y))
        tau <- (2 * S)/(n * (n - 1))
        slopes <- unlist(lapply(index, function(i, x, y) (y[i] - 
            y[1:(i - 1)])/(x[i] - x[1:(i - 1)]), x, y))
        slopes <- sort(slopes[is.finite(slopes)])
        slope <- median(slopes)
        intercept <- median(y) - slope * median(x)
        estimate <- c(tau, slope, intercept)
        names(estimate) <- c("tau", "slope", "intercept")
        method <- "Kendall's Test for Trend"
        estimation.method <- paste("slope:      Theil/Sen Estimator", 
            "intercept:  Conover's Estimator", sep = paste("\n", 
                space(33), sep = ""))
        var.S <- vark(x, y)
        if (correct) {
            method <- paste(method, "(with continuity correction)", 
                sep = paste("\n", space(33), sep = ""))
            stat <- (S - sign(S))/sqrt(var.S)
        }
        else stat <- S/sqrt(var.S)
        names(stat) <- "z"
        p.value <- switch(alternative, greater = 1 - pnorm(stat), 
            less = pnorm(stat), two.sided = 2 * pnorm(-abs(stat)))
    }
    parameters <- NULL
    null.value <- 0
    attr(null.value, "names") <- "tau"
    data.name <- c(data.name, data.name.x)
    ret.list <- list(statistic = stat, parameters = parameters, 
        p.value = p.value, estimate = estimate, null.value = null.value, 
        alternative = alternative, method = method, estimation.method = estimation.method, 
        sample.size = n, data.name = data.name, bad.obs = bad.obs, 
        S = S, var.S = var.S, slopes = slopes)
    if (!miss && ci.slope) {
        if (n < 3) {
            stop(paste("When ci.slope=TRUE, there must be at least", 
                "3 non-missing, finite observations"))
        }
        N.prime <- length(slopes)
        type <- switch(alternative, two.sided = "two-sided", 
            greater = "lower", less = "upper")
        alpha <- 1 - conf.level
        Z <- ifelse(type == "two-sided", qnorm(1 - alpha/2), 
            qnorm(conf.level))
        C.alpha <- Z * sqrt(var.S)
        M1 <- (N.prime - C.alpha)/2
        M2 <- (N.prime + C.alpha)/2
        limits <- switch(type, `two-sided` = approx(1:N.prime, 
            slopes, xout = c(M1, M2 + 1))$y, lower = c(approx(1:N.prime, 
            slopes, xout = M1)$y, Inf), upper = c(-Inf, approx(1:N.prime, 
            slopes, xout = M2 + 1)$y))
        names(limits) <- c("LCL", "UCL")
        if (any(is.na(limits))) 
            warning(paste("Sample size too small for Normal approximation", 
                "for confidence interval for slope.\n"))
        interval <- list(name = "Confidence", parameter = "slope", 
            limits = limits, type = type, method = paste("Gilbert's Modification", 
                "of Theil/Sen Method", sep = paste("\n", space(33), 
                  sep = "")), conf.level = conf.level, sample.size = N.prime)
        oldClass(interval) <- "intervalEstimate"
        ret.list <- c(ret.list, list(interval = interval))
    }
    if (!is.null(parent.of.data)) 
        ret.list$parent.of.data <- parent.of.data
    if (!is.null(subset.expression)) 
        ret.list$subset.expression <- subset.expression
    oldClass(ret.list) <- "htest"
    ret.list
}
