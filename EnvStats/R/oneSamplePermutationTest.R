oneSamplePermutationTest <-
function (x, alternative = "two.sided", mu = 0, exact = FALSE, 
    n.permutations = 5000, seed = NULL, ...) 
{
    if (!is.numeric(x)) 
        stop("'x' must be a numeric vector")
    data.name <- deparse(substitute(x))
    if ((bad.obs <- sum(!(x.ok <- is.finite(x)))) > 0) {
        is.not.finite.warning(x)
        x <- x[x.ok]
        warning(paste(bad.obs, "observations with NA/NaN/Inf in 'x' removed."))
    }
    nx <- length(x)
    if (nx < 2) 
        stop("'x' must contain at least 2 non-missing observations")
    alternative <- match.arg(alternative, c("two.sided", "less", 
        "greater"))
    if (!is.vector(mu, mode = "numeric") || is.factor(mu) || 
        length(mu) != 1 || !is.finite(mu)) 
        stop("'mu' must be a single finite numeric value")
    if (!exact) {
        if (!is.vector(n.permutations, mode = "numeric") || is.factor(n.permutations) || 
            length(n.permutations) != 1 || n.permutations != 
            trunc(n.permutations) || n.permutations < 1) 
            stop("'n.permutations' must be a positive integer")
        if (!is.null(seed) && (!is.numeric(seed) || length(seed) != 
            1 || seed != trunc(seed) || seed < 0 || seed > 1000)) 
            stop("'seed' must be an integer between 0 and 1000")
    }
    dff <- x - mu
    method <- "One-Sample Permutation Test"
    estimate <- mean(x)
    names(estimate) <- "Mean"
    if (exact) {
        if (nx <= 15) {
            sign.mat <- t(permute.signs(nx))
            stat.dist <- apply(sign.mat * abs(dff), 2, sum)
        }
        else {
            n.perms <- 2^nx
            stat.dist <- numeric(n.perms)
            sign.vec <- numeric(nx)
            for (i in 1:n.perms) {
                sign.vec[1:nx] <- base(i - 1, base = 2, num.digits = nx)
                sign.vec[sign.vec == 0] <- -1
                stat.dist[i] <- sum(sign.vec * abs(dff))
            }
        }
        method <- paste(method, "\n", space(33), "(Exact)", sep = "")
    }
    else {
        set.seed(seed)
        if ((nx * n.permutations) <= 15 * (2^15)) {
            sign.mat <- matrix(sample(c(-1, 1), size = nx * n.permutations, 
                replace = TRUE), nrow = nx, ncol = n.permutations)
            stat.dist <- apply(sign.mat * abs(dff), 2, sum)
        }
        else {
            stat.dist <- numeric(n.permutations)
            sign.vec <- numeric(nx)
            for (i in 1:n.permutations) {
                sign.vec[1:nx] <- sample(c(-1, 1), size = nx, 
                  replace = TRUE)
                stat.dist[i] <- sum(sign.vec * abs(dff))
            }
        }
        method <- paste(method, "\n", space(33), "(Based on Sampling", 
            "\n", space(33), "Permutation Distribution", "\n", 
            space(33), n.permutations, " Times)", sep = "")
    }
    if (alternative == "two.sided") {
        stat.dist <- abs(stat.dist)
        stat <- abs(sum(dff))
        names(stat) <- ifelse(mu == 0, "|Sum(x)|", paste("|Sum(x - ", 
            mu, ")|", sep = ""))
        p.value <- mean(stat.dist >= stat)
    }
    else {
        stat <- sum(dff)
        names(stat) <- ifelse(mu == 0, "Sum(x)", paste("Sum(x - ", 
            format(mu, nsmall = 0, ...), ")", sep = ""))
        if (alternative == "less") 
            p.value <- mean(stat.dist <= stat)
        else p.value <- mean(stat.dist >= stat)
    }
    parameters <- NULL
    null.value <- mu
    names(null.value) <- "Mean (Median)"
    sample.size <- c(n = nx)
    ret.list <- list(statistic = stat, parameters = parameters, 
        p.value = p.value, estimate = estimate, null.value = null.value, 
        alternative = alternative, method = method, estimation.method = NULL, 
        sample.size = sample.size, data.name = data.name, bad.obs = bad.obs, 
        stat.dist = stat.dist, exact = exact)
    if (!exact) 
        ret.list <- c(ret.list, list(seed = seed))
    oldClass(ret.list) <- "permutationTest"
    ret.list
}
