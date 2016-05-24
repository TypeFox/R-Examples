eqgammaAlt <-
function (x, p = 0.5, method = "mle", ci = FALSE, ci.type = "two-sided", 
    conf.level = 0.95, normal.approx.transform = "kulkarni.powar", 
    digits = 0) 
{
    if (!is.vector(p, mode = "numeric")) 
        stop("'p' must be a numeric vector.")
    if (any(!is.finite(p))) 
        stop("NA/NaN/Inf values not allowed in 'p'.")
    if (any(p < 0) || any(p > 1)) 
        stop("All values of 'p' must be between 0 and 1")
    method <- match.arg(method, c("mle", "bcmle", "mme", "mmue"))
    x.is.est.obj <- data.class(x) == "estimate" || data.class(x) == 
        "estimateCensored"
    if (x.is.est.obj) {
        if (ci) 
            stop(paste("When ci=TRUE the argument 'x' must be the raw data,", 
                "not the result of parameter estimation."))
        if (x$distribution != "Gamma") 
            stop(paste("'eqgamma' estimates quantiles for a gamma distribution. ", 
                "You have supplied an object that assumes a different distribution."))
        class.x <- oldClass(x)
        if (!is.null(x$interval)) {
            x <- x[-match("interval", names(x))]
            oldClass(x) <- class.x
        }
        parameters <- x$parameters
        if (names(parameters)[1] == "shape") 
            stop(paste("You supplied the result of parameter estimation", 
                "for the usual parameterization of the Gamma distribution. ", 
                "Use the function eqgamma."))
        shape <- parameters["cv"]^-2
        scale <- parameters["mean"]/shape
        n <- x$sample.size
        ret.list <- x
    }
    else {
        if (!is.vector(x, mode = "numeric")) 
            stop(paste("'x' must be either a list that inherits from", 
                "the class 'estimate', or else a numeric vector"))
        data.name <- deparse(substitute(x))
        if ((bad.obs <- sum(!(x.ok <- is.finite(x)))) > 0) {
            is.not.finite.warning(x)
            x <- x[x.ok]
            warning(paste(bad.obs, "observations with NA/NaN/Inf in 'x' removed."))
        }
        n <- length(x)
        if (n < 2 || any(x < 0) || length(unique(x)) < 2) 
            stop(paste("'x' must contain at least 2 non-missing distinct values,", 
                "and all non-missing values of x must be non-negative. ", 
                "This is not true for 'x' =", data.name))
        ret.list <- egammaAlt(x, method = method)
        ret.list$data.name <- data.name
        ret.list$bad.obs <- bad.obs
        parameters <- ret.list$parameters
        shape <- parameters["cv"]^-2
        scale <- parameters["mean"]/shape
    }
    q <- qgamma(p, shape = shape, scale = scale)
    if (length(p) == 1 && p == 0.5) 
        names(q) <- "Median"
    else {
        pct <- round(100 * p, digits)
        names(q) <- paste(pct, number.suffix(pct), " %ile", sep = "")
    }
    ret.list <- c(ret.list, list(quantiles = q))
    ret.list$quantile.method <- paste("Quantile(s) Based on\n", 
        space(33), ret.list$method, sep = "")
    if (x.is.est.obj) 
        oldClass(ret.list) <- class.x
    else oldClass(ret.list) <- "estimate"
    if (ci) {
        if (length(p) > 1 || p <= 0 || p >= 1) 
            stop(paste("When 'ci' = TRUE, 'p' must be a scalar", 
                "larger than 0 and less than 1."))
        ci.type <- match.arg(ci.type, c("two-sided", "lower", 
            "upper"))
        normal.approx.transform <- match.arg(normal.approx.transform, 
            c("kulkarni.powar", "cube.root", "fourth.root"))
        switch(normal.approx.transform, kulkarni.powar = {
            p.trans <- ifelse(shape > 1.5, 0.246, -0.0705 - 0.178 * 
                shape + 0.475 * sqrt(shape))
            string <- paste("Kulkarni & Powar (2010)\n", space(33), 
                "transformation to Normality\n", space(33), "based on ", 
                method, " of 'shape'", sep = "")
        }, cube.root = {
            p.trans <- 1/3
            string <- paste("Wilson & Hilferty (1931) cube-root\n", 
                space(33), "transformation to Normality", sep = "")
        }, fourth.root = {
            p.trans <- 1/4
            string <- paste("Hawkins & Wixley (1986) fourth-root\n", 
                space(33), "transformation to Normality", sep = "")
        })
        Y <- x^p.trans
        ci.ret.list <- eqnorm(Y, p = p, ci = TRUE, ci.type = ci.type, 
            conf.level = conf.level, digits = digits)$interval
        limits <- ci.ret.list$limits
        if (ci.type == "upper") 
            limits["LCL"] <- 0
        if (ci.type %in% c("two-sided", "lower") & limits["LCL"] < 
            0) {
            limits["LCL"] <- 0
            warning(paste("Normal approximation to Gamma distribution", 
                "to compute lower confidence limit of quantile", 
                "not accurate for this case"))
        }
        ci.ret.list$limits <- limits^(1/p.trans)
        ci.ret.list$method <- paste(ci.ret.list$method, " using\n", 
            space(33), string, sep = "")
        ret.list$interval <- ci.ret.list
    }
    ret.list
}
