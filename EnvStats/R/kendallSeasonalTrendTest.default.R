kendallSeasonalTrendTest.default <-
function (y, season, year, alternative = "two.sided", correct = TRUE, 
    ci.slope = TRUE, conf.level = 0.95, independent.obs = TRUE, 
    data.name = NULL, season.name = NULL, year.name = NULL, parent.of.data = NULL, 
    subset.expression = NULL, ...) 
{
    if (is.null(data.name)) 
        data.name <- deparse(substitute(y))
    y <- as.vector(unlist(y))
    if (!is.numeric(y)) 
        stop("All elements of 'y' must be numeric")
    if (is.null(season.name)) 
        season.name <- deparse(substitute(season))
    data.name <- c(data.name, season.name, year.name)
    names(data.name) <- c("y", "season", "year")[1:length(data.name)]
    n <- length(y)
    if (missing(season) || missing(year)) 
        stop("When 'y' is a vector you must supply both 'season' and 'year'")
    if (length(season) != n || length(year) != n) 
        stop("'season' and 'year' must both be the same length as 'y'")
    if (!(is.numeric(season) || is.factor(season) || is.character(season))) 
        stop("'season' must be a numeric or character vector or a factor")
    if (is.numeric(season) && !all(is.finite(season))) 
        stop(paste("Missing (NA), infinite (-Inf, Inf),", "and undefined (Nan) values not allowed in", 
            "'season'"))
    if (!is.numeric(year)) 
        stop("'year' must be a numeric vector")
    if (!all(is.finite(year))) 
        stop(paste("Missing (NA), infinite (-Inf, Inf),", "and undefined (Nan) values not allowed in", 
            "'year'"))
    if (!is.numeric(season)) {
        season <- as.character(season)
        season.names <- unique(season)
        season <- match(season, season.names)
    }
    else season.names <- as.character(sort(unique(season)))
    if (!independent.obs && !all(table(season, year) == 1)) 
        stop(paste("When independent.obs=FALSE,", "there must be one observation (possibly NA)", 
            "per season per year"))
    finite.index <- is.finite(y)
    y.no.na <- y[finite.index]
    season.no.na <- season[finite.index]
    year.no.na <- year[finite.index]
    n.yrs <- n.unique(year.no.na)
    if (n.yrs < 2) 
        stop("There must be at least 2 years of non-missing data.")
    n.seasons <- n.unique(season.no.na)
    if (n.seasons == 1) 
        stop("Only one season.  Use the function 'kendallTrendTest'.")
    year.by.season <- split(year.no.na, season.no.na)
    y.by.season <- split(y.no.na, season.no.na)
    season.names <- season.names[as.numeric(names(y.by.season))]
    n.vec <- sapply(y.by.season, length)
    if (all(sapply(year.by.season, n.unique) < 2)) 
        stop(paste("There must be observations in at least", 
            "2 separate years for at least one season."))
    alternative <- match.arg(alternative, c("two.sided", "greater", 
        "less"))
    if (ci.slope && !is.vector(conf.level, mode = "numeric") || 
        length(conf.level) != 1 || conf.level <= 0 || conf.level >= 
        1) 
        stop("'conf.level' must be a numeric scalar between 0 and 1")
    estimate.mat <- matrix(0, n.seasons, 3)
    S.vec <- numeric(n.seasons)
    var.S.vec <- numeric(n.seasons)
    slopes <- numeric(0)
    for (i in 1:n.seasons) {
        S.list <- kendallTrendTest(y.by.season[[i]], year.by.season[[i]], 
            ci.slope = FALSE, warn = FALSE)
        estimate.mat[i, ] <- S.list$estimate
        S.vec[i] <- S.list$S
        var.S.vec[i] <- S.list$var.S
        slopes <- c(slopes, S.list$slopes)
    }
    names(S.vec) <- season.names
    names(var.S.vec) <- season.names
    dimnames(estimate.mat) <- list(season.names, c("tau", "slope", 
        "intercept"))
    slopes <- sort(slopes, na.last = NA)
    na.index <- is.na(S.vec)
    n.seasons.actual <- n.seasons - sum(na.index)
    tau <- sum(n.vec[!na.index] * estimate.mat[!na.index, "tau"])/sum(n.vec[!na.index])
    slope <- median(slopes)
    intercept <- median(estimate.mat[!na.index, "intercept"])
    estimate <- c(tau, slope, intercept)
    names(estimate) <- c("tau", "slope", "intercept")
    method <- "Seasonal Kendall Test for Trend"
    sep.string <- paste("\n", space(33), sep = "")
    estimation.method <- paste("tau:        Weighted Average of", 
        "            Seasonal Estimates", "slope:      Hirsch et al.'s", 
        "            Modification of", "            Thiel/Sen Estimator", 
        "intercept:  Median of", "            Seasonal Estimates", 
        sep = sep.string)
    var.cov.S <- diag(var.S.vec)
    dimnames(var.cov.S) <- list(season.names, season.names)
    if (independent.obs) {
        z.vec <- S.vec[!na.index]/sqrt(var.S.vec[!na.index])
        chisq <- sum(z.vec^2) - n.seasons.actual * mean(z.vec)^2
        p.chisq <- 1 - pchisq(chisq, df = n.seasons.actual - 
            1)
    }
    else {
        method <- paste(method, "Modified for Serial Correlation", 
            sep = sep.string)
        y.by.season <- split(y, season)
        rank.mat <- sapply(y.by.season, rank.w.na, na.action = "mid.rank", 
            warn = FALSE)
        K.vec <- numeric((n.seasons * (n.seasons - 1))/2)
        R.vec <- K.vec
        count <- 1
        for (i in 2:n.seasons) {
            for (j in 1:(i - 1)) {
                K.vec[count] <- kendallTrendTest(rank.mat[, i], 
                  rank.mat[, j], ci.slope = FALSE, warn = FALSE)$S
                R.vec[count] <- 4 * sum(rank.mat[, i] * rank.mat[, 
                  j]) - n.yrs * (n.vec[i] + 1) * (n.vec[j] + 
                  1)
                count <- count + 1
            }
        }
        cov.S.vec <- (K.vec + R.vec)/3
        count <- 1
        for (i in 2:n.seasons) {
            for (j in 1:(i - 1)) {
                var.cov.S[i, j] <- cov.S.vec[count]
                var.cov.S[j, i] <- var.cov.S[i, j]
                count <- count + 1
            }
        }
        C.mat <- cbind(1, -diag(n.seasons - 1))
        m <- diag(2/(n.vec * (n.vec - 1)))
        tau.vec <- estimate.mat[, "tau"]
        vec <- C.mat %*% tau.vec
        var.vec <- C.mat %*% (m %*% var.cov.S %*% t(m)) %*% t(C.mat)
        if (qr(var.vec)$rank < (n.seasons - 1)) {
            chisq <- NA
            p.chisq <- NA
            warning(paste("Could not perform Pseudo-Heterogeneity Test", 
                "because of singluarity in variance matrix"))
        }
        else {
            chisq <- t(vec) %*% solve(var.vec) %*% vec
            p.chisq <- 1 - pchisq(chisq, df = n.seasons - 1)
        }
    }
    S <- sum(S.vec, na.rm = TRUE)
    var.S <- sum(var.cov.S, na.rm = TRUE)
    if (correct) {
        method <- paste(method, sep.string, "(with continuity correction)", 
            sep = "")
        z <- (S - sign(S))/sqrt(var.S)
    }
    else z <- S/sqrt(var.S)
    p.z <- switch(alternative, greater = 1 - pnorm(z), less = pnorm(z), 
        two.sided = 2 * pnorm(-abs(z)))
    stat <- c(chisq, z)
    names(stat) <- c("Chi-Square (Het)", "z (Trend)")
    p.value <- c(p.chisq, p.z)
    names(p.value) <- names(stat)
    parameters <- c(df = n.seasons.actual - 1)
    null.value <- rep(0, n.seasons)
    names(null.value) <- rep("tau", n.seasons)
    alternative.chisq <- paste("The seasonal taus are not all equal", 
        "(Chi-Square Heterogeneity Test)", sep = sep.string)
    alternative.z <- switch(alternative, greater = paste("At least one seasonal tau > 0", 
        "(z Trend Test)", sep = sep.string), less = paste("At least one seasonal tau < 0", 
        "(z Trend Test)", sep = sep.string), two.sided = paste("At least one seasonal tau != 0", 
        "and all non-zero tau's have the", "same sign (z Trend Test)", 
        sep = sep.string))
    alt <- paste(alternative.chisq, alternative.z, sep = sep.string)
    n.vec <- c(n.vec, Total = sum(n.vec))
    names(n.vec) <- c(season.names, "Total")
    ret.list <- list(statistic = stat, parameters = parameters, 
        p.value = p.value, estimate = estimate, null.value = null.value, 
        alternative = alt, method = method, estimation.method = estimation.method, 
        sample.size = n.vec, data.name = data.name, bad.obs = sum(!finite.index), 
        seasonal.S = S.vec)
    if (independent.obs) 
        ret.list <- c(ret.list, list(var.seasonal.S = var.S.vec, 
            seasonal.estimates = estimate.mat))
    else ret.list <- c(ret.list, list(var.cov.seasonal.S = var.cov.S, 
        seasonal.estimates = estimate.mat))
    if (ci.slope) {
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
            limits = limits, type = type, method = paste("Gilbert's Modification of", 
                "Theil/Sen Method", sep = sep.string), conf.level = conf.level, 
            sample.size = N.prime)
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
