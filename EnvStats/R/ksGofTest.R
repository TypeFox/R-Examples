ksGofTest <-
function (x, y = NULL, alternative = "two.sided", exact = NULL, 
    distribution = "norm", param.list = NULL, estimate.params = ifelse(is.null(param.list), 
        TRUE, FALSE), est.arg.list = NULL, n.param.est = NULL, 
    digits = .Options$digits, data.name.x = NULL, data.name.y = NULL) 
{
    alt.expanded <- if (!missing(alternative)) 
        char.expand(alternative, c("two.sided", "greater", "less"), 
            stop("argument 'alternative' must match one of \"greater\",\n                                                    \"less\", \"two.sided\"."))
    else alternative
    if (!is.vector(x, mode = "numeric") || is.factor(x)) 
        stop("'x' must be a numeric vector")
    if (is.null(data.name.x)) 
        data.name.x <- deparse(substitute(x))
    if (is.null(y)) {
        if ((bad.obs <- sum(!(x.ok <- is.finite(x)))) > 0) {
            is.not.finite.warning(x)
            x <- x[x.ok]
            warning(paste(bad.obs, "observations with NA/NaN/Inf in 'x' removed."))
        }
        names(x) <- NULL
        if (estimate.params) 
            check.da.list <- check.distribution.args(distribution, 
                check.params = FALSE)
        else {
            if (is.null(param.list)) 
                stop(paste("When 'estimate.params=F' you must supply", 
                  "the argument 'param.list'"))
            check.da.list <- check.distribution.args(distribution, 
                param.list)
        }
        dist.abb <- check.da.list$dist.abb
        dist.name <- check.da.list$dist.name
        dist.type <- check.da.list$dist.type
        if (dist.type != "Continuous") 
            stop(paste("The Kolmogorov-Smirnov goodness-of-fit test", 
                "is currenlty implemented only for continuous distributions"))
        n.dist.params <- check.da.list$n.dist.params
        dist.params.names <- check.da.list$dist.params.names
        if (estimate.params) {
            if (EnvStats::Distribution.df[dist.abb, "Estimation.Method(s)"] == 
                "") 
                stop(paste("No estimation method for", dist.name))
            warning(paste("The standard Kolmogorov-Smirnov test is very", 
                "conservative", "(Type I error smaller than assumed; high Type II error) for", 
                "testing departures from the", dist.name, "distribution when", 
                "you have to estimate the distribution parameters.\n\n"))
            ename <- paste("e", dist.abb, sep = "")
            est.list <- do.call(ename, c(list(x = x), est.arg.list))
            dist.params <- est.list$parameters
            estimation.method <- est.list$method
            param.list <- as.list(dist.params)
            names.param.list <- names(param.list) <- names(dist.params)
        }
        else {
            param.list <- check.da.list$param.list
            dist.params <- unlist(param.list)
            estimation.method <- NULL
        }
        if (is.null(n.param.est)) 
            n.param.est <- ifelse(estimate.params, n.dist.params, 
                0)
        nx <- length(x)
        pname <- paste("p", dist.abb, sep = "")
        test.list <- do.call("ks.test", c(list(x = x, y = pname, 
            alternative = alternative, exact = exact), param.list))
        if (n.param.est == 0) {
            if (any(dist.abb == c("beta", "chisq", "f")) && param.list$ncp == 
                0) 
                hyp.dist <- paste(dist.name, "(", paste(paste(dist.params.names[-n.dist.params], 
                  signif(dist.params[-n.dist.params], digits = digits), 
                  sep = " = "), collapse = ", "), ")", sep = "")
            else hyp.dist <- paste(dist.name, "(", paste(paste(dist.params.names, 
                signif(dist.params, digits = digits), sep = " = "), 
                collapse = ", "), ")", sep = "")
            alt.hyp <- paste("True cdf ", ifelse(alt.expanded == 
                "two.sided", "does not equal ", paste("is", alt.expanded, 
                "than ")), "the\n", space(33), hyp.dist, "\n", 
                space(33), "Distribution.", sep = "")
        }
        else {
            hyp.dist <- dist.name
            alt.hyp <- paste("True cdf ", ifelse(alt.expanded == 
                "two.sided", "does not equal ", paste("is", alt.expanded, 
                "than ")), "the\n", space(33), hyp.dist, " Distribution.", 
                sep = "")
        }
        ret.val <- list(distribution = hyp.dist, dist.abb = dist.abb, 
            distribution.parameters = dist.params, n.param.est = n.param.est, 
            estimation.method = estimation.method, statistic = test.list$statistic, 
            sample.size = nx, parameters = c(n = nx), p.value = test.list$p.value, 
            alternative = alt.hyp, method = "Kolmogorov-Smirnov GOF", 
            data = x, data.name = data.name.x)
        if (bad.obs > 0) 
            ret.val$bad.obs <- bad.obs
        oldClass(ret.val) <- "gof"
    }
    else {
        if (!is.vector(y, mode = "numeric") || is.factor(y)) 
            stop("'y' must be a numeric vector")
        if (is.null(data.name.y)) 
            data.name.y <- deparse(substitute(y))
        if ((alt.expanded == "less") || (alt.expanded == "more")) 
            stop("Only the two sided alternative is tested for two samples.")
        if ((bad.obs.x <- sum(!(x.ok <- is.finite(x)))) > 0) {
            is.not.finite.warning(x)
            x <- x[x.ok]
            warning(paste(bad.obs.x, "observations with NA/NaN/Inf in 'x' removed."))
        }
        names(x) <- NULL
        nx <- length(x)
        if ((bad.obs.y <- sum(!(y.ok <- is.finite(y)))) > 0) {
            is.not.finite.warning(y)
            y <- y[y.ok]
            warning(paste(bad.obs.y, "observations with NA/NaN/Inf in 'y' removed."))
        }
        names(y) <- NULL
        ny <- length(y)
        test.list <- ks.test(x, y, alternative = alternative, 
            exact = exact)
        data.list <- list(x, y)
        names(data.list) <- c(data.name.x, data.name.y)
        ret.val <- list(distribution = "Equal", statistic = test.list$statistic, 
            sample.size = c(n.x = nx, n.y = ny), parameters = c(n = nx, 
                m = ny), p.value = test.list$p.value, alternative = paste("The cdf of '", 
                data.name.x, "' does not equal\n", space(33), 
                "the cdf of '", data.name.y, "'.", sep = ""), 
            method = "2-Sample K-S GOF", data = data.list, data.name = paste(data.name.x, 
                " and ", data.name.y, sep = ""))
        if (bad.obs.x > 0 | bad.obs.y > 0) 
            ret.val$bad.obs <- paste(bad.obs.x, "and", bad.obs.y)
        oldClass(ret.val) <- "gofTwoSample"
    }
    names(ret.val$statistic) <- "ks"
    return(ret.val)
}
