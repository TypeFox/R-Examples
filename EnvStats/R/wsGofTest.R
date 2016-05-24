wsGofTest <-
function (x, method = "normal scores", alternative = "greater") 
{
    if (!is.vector(x, mode = "numeric")) 
        stop("'x' must be a numeric vector")
    data.name <- deparse(substitute(x))
    if ((bad.obs <- sum(!(x.ok <- is.finite(x)))) > 0) {
        is.not.finite.warning(x)
        x <- x[x.ok]
        warning(paste(bad.obs, "observations with NA/NaN/Inf in 'x' removed."))
    }
    n <- length(x)
    if (any(x < 0) || any(x > 1) || n < 2) 
        stop(paste("'x' must be a numeric vector with", "at least 2 non-missing elements,", 
            "and all values between 0 and 1"))
    sep.string <- paste("\n", space(33), sep = "")
    method <- match.arg(method, c("normal scores", "chi-square scores"))
    alternative <- match.arg(alternative, c("less", "greater"))
    dist <- "Uniform [0, 1]"
    dist.params <- c(min = 0, max = 1)
    if (method == "normal scores") {
        scores <- qnorm(x)
        stat <- sum(scores)/sqrt(n)
        names(stat) <- "z (G)"
        p.value <- pnorm(stat)
        names(p.value) <- names(stat)
        parameters <- NULL
    }
    else {
        scores <- -2 * log(x)
        stat <- sum(scores)
        names(stat) <- "Chi-square (C)"
        df <- 2 * n
        p.value <- pchisq(stat, df = df, lower.tail = FALSE)
        names(p.value) <- names(stat)
        parameters <- c(df = df)
    }
    if (alternative == "less") 
        p.value <- 1 - p.value
    alt <- paste("True cdf is ", alternative, " than the", sep.string, 
        dist, sep.string, "Distribution.", sep = "")
    method <- ifelse(method == "normal scores", "Normal Scores", 
        "Chi-Square Scores")
    ret.list <- list(distribution = dist, dist.abb = "unif", 
        distribution.parameters = dist.params, n.param.est = 0, 
        estimation.method = NULL, statistic = stat, sample.size = n, 
        parameters = parameters, p.value = p.value, alternative = alt, 
        method = paste("Wilk-Shapiro GOF (", method, ")", sep = ""), 
        data = x, data.name = data.name, bad.obs = bad.obs, scores = scores)
    oldClass(ret.list) <- "gof"
    ret.list
}
