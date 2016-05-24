predIntGammaAltSimultaneous <-
function (x, n.transmean = 1, k = 1, m = 2, r = 1, rule = "k.of.m", 
    delta.over.sigma = 0, pi.type = "upper", conf.level = 0.95, 
    K.tol = 1e-07, est.method = "mle", normal.approx.transform = "kulkarni.powar") 
{
    rule <- match.arg(rule, c("k.of.m", "CA", "Modified.CA"))
    pi.type <- match.arg(pi.type, c("upper", "lower"))
    est.method <- match.arg(est.method, c("mle", "bcmle", "mme", 
        "mmue"))
    normal.approx.transform <- match.arg(normal.approx.transform, 
        c("kulkarni.powar", "cube.root", "fourth.root"))
    switch(rule, k.of.m = {
        if (!is.vector(n.transmean, mode = "numeric") || length(n.transmean) != 
            1 || n.transmean != trunc(n.transmean) || n.transmean < 
            1 || !is.vector(k, mode = "numeric") || length(k) != 
            1 || k != trunc(k) || k < 1 || !is.vector(m, mode = "numeric") || 
            length(m) != 1 || m != trunc(m) || m < 1 || !is.vector(r, 
            mode = "numeric") || length(r) != 1 || r != trunc(r) || 
            r < 1 || k > m) stop(paste("'n.transmean', 'k', 'm', and 'r' must be positive integers,", 
            "and 'k' must be between 1 and 'm'"))
    }, CA = {
        if (!is.vector(n.transmean, mode = "numeric") || length(n.transmean) != 
            1 || n.transmean != trunc(n.transmean) || n.transmean < 
            1 || !is.vector(m, mode = "numeric") || length(m) != 
            1 || m != trunc(m) || m < 1 || !is.vector(r, mode = "numeric") || 
            length(r) != 1 || r != trunc(r) || r < 1) stop("'n.transmean', 'm', and 'r' must be positive integers")
    }, Modified.CA = {
        if (!is.vector(n.transmean, mode = "numeric") || length(n.transmean) != 
            1 || n.transmean != trunc(n.transmean) || n.transmean < 
            1 || !is.vector(m, mode = "numeric") || length(m) != 
            1 || m != trunc(m) || m < 1 || !is.vector(r, mode = "numeric") || 
            length(r) != 1 || r != trunc(r) || r < 1) stop("'n.transmean', 'm', and 'r' must be positive integers")
        m <- 4
    })
    if (!is.vector(delta.over.sigma, mode = "numeric") || length(delta.over.sigma) != 
        1 || !is.finite(delta.over.sigma)) 
        stop("'delta.over.sigma' must be a finite numeric scalar.")
    if (!is.vector(conf.level, mode = "numeric") || length(conf.level) != 
        1 || conf.level <= 0 || conf.level >= 1) 
        stop("'conf.level' must be a scalar greater than 0 and less than 1.")
    if (!is.vector(x, mode = "numeric")) 
        stop("'x' must be a numeric vector")
    data.name <- deparse(substitute(x))
    if ((bad.obs <- sum(!(x.ok <- is.finite(x)))) > 0) {
        is.not.finite.warning(x)
        x <- x[x.ok]
        warning(paste(bad.obs, "observations with NA/NaN/Inf in 'x' removed."))
    }
    if (any(x < 0)) 
        stop("All non-missing values of 'x' must be non-negative")
    n <- length(x)
    if (n < 2 || length(unique(x)) < 2) 
        stop(paste("'x' must contain at least 2 non-missing distinct values. ", 
            "This is not true for 'x' =", data.name))
    dum.list <- egamma(x = x, method = est.method)
    shape <- dum.list$parameters["shape"]
    scale <- dum.list$parameters["scale"]
    switch(normal.approx.transform, kulkarni.powar = {
        p <- ifelse(shape > 1.5, 0.246, -0.0705 - 0.178 * shape + 
            0.475 * sqrt(shape))
        string <- paste("Kulkarni & Powar (2010)\n", space(33), 
            "transformation to Normality\n", space(33), "based on ", 
            dum.list$method, " of 'shape'", sep = "")
    }, cube.root = {
        p <- 1/3
        string <- paste("Wilson & Hilferty (1931) cube-root\n", 
            space(33), "transformation to Normality", sep = "")
    }, fourth.root = {
        p <- 1/4
        string <- paste("Hawkins & Wixley (1986) fourth-root\n", 
            space(33), "transformation to Normality", sep = "")
    })
    names(p) <- NULL
    Y <- x^p
    ret.list <- predIntNormSimultaneous(Y, n.mean = n.transmean, 
        k = k, m = m, r = r, rule = rule, delta.over.sigma = delta.over.sigma, 
        pi.type = pi.type, conf.level = conf.level, K.tol = K.tol)
    ret.list$data.name <- data.name
    ret.list$bad.obs <- bad.obs
    dum.list <- egammaAlt(x = x, method = est.method)
    ret.list$parameters <- dum.list$parameters
    ret.list$method <- dum.list$method
    ret.list$distribution <- "Gamma"
    names(ret.list$interval)[names(ret.list$interval) == "n.mean"] <- "n.transmean"
    ret.list$interval$method <- paste(ret.list$interval$method, 
        " using\n", space(33), string, sep = "")
    limits <- ret.list$interval$limits
    if (pi.type == "upper") 
        limits["LPL"] <- 0
    if (pi.type %in% c("two-sided", "lower") & limits["LPL"] < 
        0) {
        limits["LPL"] <- 0
        warning("Normal approximation not accurate for this case")
    }
    ret.list$interval$limits <- limits^(1/p)
    ret.list$interval$normal.transform.power <- p
    ret.list
}
