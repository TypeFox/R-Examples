chisqGofTest <-
function (x, n.classes = NULL, cut.points = NULL, distribution = "norm", 
    param.list = NULL, estimate.params = ifelse(is.null(param.list), 
        TRUE, FALSE), est.arg.list = NULL, n.param.est = NULL, 
    correct = NULL, digits = .Options$digits) 
{
    if (!is.vector(x, mode = "numeric") || is.factor(x)) 
        stop("'x' must be a numeric vector")
    data.name <- deparse(substitute(x))
    if ((bad.obs <- sum(!(x.ok <- is.finite(x)))) > 0) {
        is.not.finite.warning(x)
        x <- x[x.ok]
        warning(paste(bad.obs, "observations with NA/NaN/Inf in 'x' removed."))
    }
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
    if (dist.name == "Stable") 
        stop("No method for Stable Distribution.")
    n.dist.params <- check.da.list$n.dist.params
    dist.params.names <- check.da.list$dist.params.names
    if (estimate.params) {
        if (EnvStats::Distribution.df[dist.abb, "Estimation.Method(s)"] == 
            "") 
            stop(paste("No estimation method for", dist.name, 
                "Distribution"))
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
    if (!is.null(cut.points)) {
        n.cut.points <- length(cut.points)
        n.classes <- n.cut.points - 1
    }
    else {
        if (is.null(n.classes)) 
            n.classes <- ceiling(2 * (length(x)^(2/5)))
    }
    if (is.null(correct)) 
        correct <- n.classes == 2
    n <- length(x)
    pdist.abb <- paste("p", dist.abb, sep = "")
    qdist.abb <- paste("q", dist.abb, sep = "")
    if (is.null(cut.points)) {
        if (dist.type != "Continuous") 
            stop("You must supply the cutpoints for a discrete distribution")
        else {
            num <- pmin(floor(1 + n.classes * do.call(pdist.abb, 
                c(list(x), param.list))), n.classes)
            prob <- rep(1/n.classes, n.classes)
            cut.points <- do.call(qdist.abb, c(list(seq(0, 1, 
                by = 1/n.classes)), param.list))
        }
    }
    else {
        n.cut.pts <- cut.points
        if (cut.points[1] == -Inf) 
            n.cut.pts[1] <- min(x, cut.points[-1]) - 1
        if (cut.points[n.cut.points] == Inf) 
            n.cut.pts[n.cut.points] <- max(x, cut.points[-n.cut.points]) + 
                1
        if (sum(!(is.finite(n.cut.pts))) > 0) 
            stop("Please remove Inf's or NA's from cut.points\n")
        num <- cut(x, n.cut.pts)
        if ((bad.obs <- sum(!(num.ok <- is.finite(num)))) > 0) {
            is.not.finite.warning(x)
            num <- num[num.ok]
            warning(paste(bad.obs, "observations do not fall within within the given cutpoints.  \n            There are removed."))
            if (length(x) - bad.obs < 2) 
                stop("Less than 2 non-missing values. Impossible to continue\n")
        }
        prob <- diff(do.call(pdist.abb, c(list(cut.points), param.list)))
    }
    count <- tabulate(num, n.classes)
    xpec <- n * prob
    X2.components <- (abs(count - xpec) - (if (correct) 
        0.5
    else 0))^2/xpec
    X2 <- sum(X2.components)
    if (any(xpec < 5)) 
        warning("Expected counts < 5. Chi-squared approximation may not\n                        be appropriate.")
    parameters <- n.classes - n.param.est - 1
    p.value <- 1 - pchisq(X2, parameters)
    if (n.param.est == 0) {
        if (any(dist.abb == c("beta", "chisq", "f")) && param.list$ncp == 
            0) 
            hyp.dist <- paste(dist.name, "(", paste(paste(dist.params.names[-n.dist.params], 
                signif(dist.params[-n.dist.params], digits = digits), 
                sep = " = "), collapse = ", "), ")", sep = "")
        else hyp.dist <- paste(dist.name, "(", paste(paste(dist.params.names, 
            signif(dist.params, digits = digits), sep = " = "), 
            collapse = ", "), ")", sep = "")
        alt.hyp <- paste("True cdf does not equal the\n", space(33), 
            hyp.dist, "\n", space(33), "Distribution.", sep = "")
    }
    else {
        hyp.dist <- dist.name
        alt.hyp <- paste("True cdf does not equal the\n", space(33), 
            hyp.dist, " Distribution.", sep = "")
    }
    ret.val <- list(distribution = hyp.dist, dist.abb = dist.abb, 
        distribution.parameters = dist.params, n.param.est = n.param.est, 
        estimation.method = estimation.method, statistic = X2, 
        sample.size = n, parameters = parameters, p.value = p.value, 
        alternative = alt.hyp, method = "Chi-square GOF", data = x, 
        data.name = data.name, bad.obs = bad.obs, cut.points = cut.points, 
        counts = count, expected = xpec, X2.components = X2.components)
    names(ret.val$statistic) <- "Chi-square"
    names(ret.val$parameters) <- "df"
    names(ret.val$p.value) <- NULL
    oldClass(ret.val) <- "gof"
    ret.val
}
