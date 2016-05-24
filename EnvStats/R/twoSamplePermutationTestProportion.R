twoSamplePermutationTestProportion <-
function (x, y, x.and.y = "Binomial Outcomes", alternative = "two.sided", 
    tol = sqrt(.Machine$double.eps)) 
{
    x.and.y <- match.arg(x.and.y, c("Binomial Outcomes", "Number Successes and Trials"))
    alternative <- match.arg(alternative, c("two.sided", "less", 
        "greater"))
    x.name <- deparse(substitute(x))
    y.name <- deparse(substitute(y))
    data.name <- c(x.name, y.name)
    if (x.and.y == "Number Successes and Trials") {
        if (length(x) != 2 || !all(is.finite(x)) || !all(x == 
            trunc(x)) || !all(x >= 0)) 
            stop(paste("When x.and.y=\"Number Successes and Trials\",", 
                "'x' must be a numeric vector", "with two elements, and both elements must be", 
                "non-negative integers"))
        if (length(y) != 2 || !all(is.finite(y)) || !all(y == 
            trunc(y)) || !all(y >= 1) || !all(y >= x)) 
            stop(paste("When x.and.y=\"Number Successes and Trials\",", 
                "'y' must be a numeric vector", "with two elements, both elements must be", 
                "positive integers, and each element of 'y' must be", 
                "at least as large as the corresponding element", 
                "of 'x'"))
        bad.obs <- NULL
        names(data.name) <- c("x", "n")
        x1 <- x[1]
        x2 <- x[2]
        n1 <- y[1]
        n2 <- y[2]
    }
    else {
        x <- as.numeric(x)
        y <- as.numeric(y)
        if ((bad.obs <- sum(!(x.ok <- is.finite(x)))) > 0) {
            is.not.finite.warning(x)
            x <- x[x.ok]
            warning(paste(bad.obs, "observations with NA/NaN/Inf in 'x' removed."))
        }
        if ((bad.obs <- sum(!(y.ok <- is.finite(y)))) > 0) {
            is.not.finite.warning(y)
            y <- y[y.ok]
            warning(paste(bad.obs, "observations with NA/NaN/Inf in 'y' removed."))
        }
        names(data.name) <- c("x", "y")
        n1 <- length(x)
        n2 <- length(y)
        if (n1 < 2 || n2 < 2) 
            stop(paste("When x.and.y=\"Binomial Outcomes\",", 
                "'x' and 'y' must each contain", "at least two non-missing observations"))
        levels.x.y <- sort(unique(c(x, y)))
        if (length(levels.x.y) > 2) {
            stop(paste("When x.and.y=\"Binomial Outcomes\",", 
                "'x' and 'y' must take only at most two unique values"))
        }
        if (length(levels.x.y) == 1) {
            x1 <- 0
            x2 <- 0
        }
        else {
            x1 <- sum(x == levels.x.y[2])
            x2 <- sum(y == levels.x.y[2])
        }
    }
    method <- paste("Two-Sample Permutation Test\n", space(33), 
        "Based on Differences in Proportions\n", space(33), "(Fisher's Exact Test)", 
        sep = "")
    estimate <- c(x1/n1, x2/n2)
    names(estimate) <- c("p.hat.x", "p.hat.y")
    n <- n1 + n2
    n.1s <- x1 + x2
    n.0s <- n - n.1s
    x.1s <- max(0, n1 - n.0s):min(n1, n.1s)
    stat.dist <- x.1s/n1 - (n.1s - x.1s)/n2
    probs.stat.dist <- dhyper(x = x.1s, m = n.1s, n = n.0s, k = n1)
    stat <- -diff(estimate)
    if (alternative == "two.sided") {
        abs.stat.dist <- sort(unique(abs(stat.dist)))
        n.temp <- length(abs.stat.dist)
        probs.abs.stat.dist <- numeric(length(n.temp))
        for (i in 1:n.temp) {
            probs.abs.stat.dist[i] <- sum(probs.stat.dist[abs(stat.dist) == 
                abs.stat.dist[i]])
        }
        stat <- abs(stat)
        names(stat) <- "|p.hat.x - p.hat.y|"
        stat.dist <- abs.stat.dist
        probs.stat.dist <- probs.abs.stat.dist
        p.value <- sum(probs.stat.dist[stat.dist >= stat - tol])
    }
    else {
        names(stat) <- "p.hat.x - p.hat.y"
        if (alternative == "less") 
            p.value <- sum(probs.stat.dist[stat.dist <= stat + 
                tol])
        else p.value <- sum(probs.stat.dist[stat.dist >= stat - 
            tol])
    }
    parameters <- NULL
    null.value <- 0
    names(null.value) <- "p.x - p.y"
    sample.size <- c(nx = n1, ny = n2)
    ret.list <- list(statistic = stat, parameters = parameters, 
        p.value = p.value, estimate = estimate, null.value = null.value, 
        alternative = alternative, method = method, estimation.method = NULL, 
        sample.size = sample.size, data.name = data.name, bad.obs = bad.obs, 
        stat.dist = stat.dist, probs.stat.dist = probs.stat.dist, 
        exact = TRUE)
    oldClass(ret.list) <- "permutationTest"
    ret.list
}
