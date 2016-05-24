twoSamplePermutationTestLocation <-
function (x, y, fcn = "mean", alternative = "two.sided", mu1.minus.mu2 = 0, 
    paired = FALSE, exact = FALSE, n.permutations = 5000, seed = NULL, 
    tol = sqrt(.Machine$double.eps)) 
{
    alternative <- match.arg(alternative, c("two.sided", "less", 
        "greater"))
    if (!is.vector(mu1.minus.mu2, mode = "numeric") || is.factor(mu1.minus.mu2) || 
        length(mu1.minus.mu2) != 1 || !is.finite(mu1.minus.mu2)) 
        stop("'mu1.minus.mu2' must be a single finite numeric value")
    if (!exact) {
        if (!is.vector(n.permutations, mode = "numeric") || is.factor(n.permutations) || 
            length(n.permutations) != 1 || n.permutations != 
            trunc(n.permutations) || n.permutations < 1) 
            stop("'n.permutations' must be a positive integer")
        if (!is.null(seed) && (!is.numeric(seed) || length(seed) != 
            1 || seed != trunc(seed) || seed < 0 || seed > 1000)) 
            stop("'seed' must be an integer between 0 and 1000")
    }
    if (!is.numeric(x)) 
        stop("'x' must be a numeric vector")
    x.name <- deparse(substitute(x))
    if (!is.numeric(y)) 
        stop("'y' must be a numeric vector")
    y.name <- deparse(substitute(y))
    if (paired) {
        if (length(y) != length(x)) 
            stop("'x' and 'y' must be the same length when paired=T")
        if ((bad.obs <- sum(!(ok <- is.finite(x) & is.finite(y)))) > 
            0) {
            is.not.finite.warning(x)
            is.not.finite.warning(y)
            x <- x[ok]
            y <- y[ok]
            warning(paste(bad.obs, "observations with NA/NaN/Inf in 'x' and 'y' removed."))
        }
    }
    else {
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
    }
    nx <- length(x)
    ny <- length(y)
    if (nx < 2 || ny < 2) 
        stop("'x' and 'y' must contain at least two non-missing observations")
    data.name <- c(x.name, y.name)
    names(data.name) <- c("x", "y")
    if (paired) {
        ret.list <- oneSamplePermutationTest(x = x - y, alternative = alternative, 
            mu = mu1.minus.mu2, exact = exact, n.permutations = n.permutations, 
            seed = seed)
        names(ret.list$estimate) <- "Mean (Median) of Differences"
        method <- "Paired-Sample Permutation Test"
        if (exact) {
            ret.list$method <- paste(method, "\n", space(33), 
                "(Exact)", sep = "")
        }
        else {
            ret.list$method <- paste(method, "\n", space(33), 
                "(Based on Sampling", "\n", space(33), "Permutation Distribution", 
                "\n", space(33), n.permutations, " Times)", sep = "")
        }
        if (ret.list$alternative == "two.sided") {
            names(ret.list$statistic) <- ifelse(mu1.minus.mu2 == 
                0, "|Sum(x-y)|", paste("|Sum(x-y) - ", mu1.minus.mu2, 
                "|", sep = ""))
        }
        else {
            names(ret.list$statistic) <- ifelse(mu1.minus.mu2 == 
                0, "Sum(x-y)", paste("Sum(x-y) -", mu1.minus.mu2))
        }
        names(ret.list$null.value) <- "Mean (Median) of Differences"
        ret.list$data.name <- data.name
    }
    else {
        fcn <- match.arg(fcn, c("mean", "median"))
        fcn.name <- fcn
        if (fcn.name == "mean") 
            fcn <- get("mean", pos = "package:base")
        else fcn <- get("median", pos = "package:stats")
        string <- ifelse(fcn.name == "mean", "Means", "Medians")
        method <- paste("Two-Sample Permutation Test\n", space(33), 
            "Based on Differences in ", string, sep = "")
        estimate <- c(fcn(x), fcn(y))
        names(estimate) <- c(paste(fcn.name, "of x"), paste(fcn.name, 
            "of y"))
        x <- x - mu1.minus.mu2
        xy <- c(x, y)
        n <- nx + ny
        if (exact) {
            if (mu1.minus.mu2 == 0 && fcn.name == "mean") {
                if (n < 20) {
                  comb.mat <- t(combn(n, nx))
                  x.comb.mat <- matrix(xy[comb.mat], nrow = choose(n, 
                    nx))
                  sum.x <- apply(x.comb.mat, 1, sum)
                  stat.dist <- sum.x * (1/nx + 1/ny) - sum(xy)/ny
                }
                else {
                  stop(paste("When fcn=\"mean\" and mu1.minus.mu2=0,", 
                    "exact method available only for", "combined sample size size less than 20"))
                }
            }
            else {
                if (n < 10) {
                  perm.mat <- permutations(n)
                  xy.perm.mat <- matrix(xy[perm.mat], nrow = factorial(n))
                  x.fcn <- apply(xy.perm.mat[, 1:nx], 1, fcn)
                  y.fcn <- apply(xy.perm.mat[, -(1:nx)], 1, fcn)
                  stat.dist <- x.fcn - y.fcn
                }
                else {
                  stop(paste("When fcn=\"median\", or when", 
                    "fcn=\"mean\" and mu1.minus.mu2!=0,", "exact method available only for", 
                    "combined sample size size less than 10"))
                }
            }
            method <- paste(method, "\n", space(33), "(Exact)", 
                sep = "")
        }
        else {
            set.seed(seed)
            if ((n * n.permutations) <= 15 * (2^15)) {
                perm.mat <- t(sapply(1:n.permutations, function(x, 
                  n) sample(n), n = n))
                xy.perm.mat <- matrix(xy[perm.mat], nrow = n.permutations)
                x.fcn <- apply(xy.perm.mat[, 1:nx], 1, fcn)
                y.fcn <- apply(xy.perm.mat[, -(1:nx)], 1, fcn)
                stat.dist <- x.fcn - y.fcn
            }
            else {
                stat.dist <- numeric(n.permutations)
                perm.vec <- 1:n
                for (i in 1:n.permutations) {
                  perm.vec[1:n] <- sample(n)
                  stat.dist[i] <- fcn(xy[perm.vec][1:nx]) - fcn(xy[perm.vec][-(1:nx)])
                }
            }
            method <- paste(method, "\n", space(33), "(Based on Sampling", 
                "\n", space(33), "Permutation Distribution", 
                "\n", space(33), n.permutations, " Times)", sep = "")
        }
        stat <- fcn(x) - fcn(y)
        if (alternative == "two.sided") {
            stat.dist <- abs(stat.dist)
            stat <- abs(stat)
            if (mu1.minus.mu2 == 0) 
                names(stat) <- paste("|", fcn.name, ".x - ", 
                  fcn.name, ".y|", sep = "")
            else names(stat) <- paste("|", fcn.name, ".x - ", 
                fcn.name, ".y - ", mu1.minus.mu2, "|", sep = "")
            p.value <- mean(stat.dist >= stat - tol)
        }
        else {
            names(stat) <- paste(fcn.name, ".x - ", fcn.name, 
                ".y", sep = "")
            if (mu1.minus.mu2 != 0) 
                names(stat) <- paste(names(stat), " - ", mu1.minus.mu2, 
                  sep = "")
            if (alternative == "less") 
                p.value <- mean(stat.dist <= stat + tol)
            else p.value <- mean(stat.dist >= stat - tol)
        }
        parameters <- NULL
        null.value <- mu1.minus.mu2
        names(null.value) <- "mu.x-mu.y"
        sample.size <- c(nx = nx, ny = ny)
        ret.list <- list(statistic = stat, parameters = parameters, 
            p.value = p.value, estimate = estimate, null.value = null.value, 
            alternative = alternative, method = method, estimation.method = NULL, 
            sample.size = sample.size, data.name = data.name, 
            bad.obs = bad.obs, stat.dist = stat.dist, exact = exact)
        if (!exact) 
            ret.list <- c(ret.list, list(seed = seed))
        oldClass(ret.list) <- "permutationTest"
    }
    ret.list
}
