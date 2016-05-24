tolIntGammaAlt <-
function (x, coverage = 0.95, cov.type = "content", ti.type = "two-sided", 
    conf.level = 0.95, method = "exact", est.method = "mle", 
    normal.approx.transform = "kulkarni.powar") 
{
    if (any(is.na(coverage))) 
        stop("Missing values not allowed in 'coverage'.")
    if (!is.numeric(coverage) || length(coverage) > 1 || coverage <= 
        0 || coverage >= 1) 
        stop("'coverage' must be a scalar greater than 0 and less than 1.")
    cov.type <- match.arg(cov.type, c("content", "expectation"))
    ti.type <- match.arg(ti.type, c("two-sided", "lower", "upper"))
    method <- match.arg(method, c("exact", "wald.wolfowitz"))
    est.method <- match.arg(est.method, c("mle", "bcmle", "mme", 
        "mmue"))
    normal.approx.transform <- match.arg(normal.approx.transform, 
        c("kulkarni.powar", "cube.root", "fourth.root"))
    if (!is.numeric(conf.level) || length(conf.level) > 1 || 
        conf.level <= 0 || conf.level >= 1) 
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
    switch(normal.approx.transform, kulkarni.powar = {
        shape <- dum.list$parameters["shape"]
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
    ret.list <- tolIntNorm(Y, coverage = coverage, cov.type = cov.type, 
        ti.type = ti.type, conf.level = conf.level, method = method)
    ret.list$data.name <- data.name
    ret.list$bad.obs <- bad.obs
    dum.list <- egammaAlt(x = x, method = est.method)
    ret.list$parameters <- dum.list$parameters
    ret.list$method <- dum.list$method
    ret.list$distribution <- "Gamma"
    ret.list$interval$method <- paste(ret.list$interval$method, 
        " using\n", space(33), string, sep = "")
    limits <- ret.list$interval$limits
    if (ti.type == "upper") 
        limits["LTL"] <- 0
    if (ti.type %in% c("two-sided", "lower") & limits["LTL"] < 
        0) {
        limits["LTL"] <- 0
        warning("Normal approximation not accurate for this case")
    }
    ret.list$interval$limits <- limits^(1/p)
    ret.list$interval$normal.transform.power <- p
    ret.list
}
