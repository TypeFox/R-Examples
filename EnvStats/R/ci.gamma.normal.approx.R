ci.gamma.normal.approx <-
function (x, shape, shape.est.method, ci.type, normal.approx.transform, 
    conf.level) 
{
    switch(normal.approx.transform, kulkarni.powar = {
        p <- ifelse(shape > 1.5, 0.246, -0.0705 - 0.178 * shape + 
            0.475 * sqrt(shape))
        method <- paste("Optimum Power Normal Approximation\n", 
            space(33), "of Kulkarni & Powar (2010)\n", space(33), 
            "using ", shape.est.method, " of 'shape'", sep = "")
    }, cube.root = {
        p <- 1/3
        method <- paste("Normal Approximation of\n", space(33), 
            "Kulkarni & Powar (2010) using\n", space(33), "Wilson & Hilferty (1931) cube-root\n", 
            space(33), "transformation to Normality", space(33), 
            "and ", shape.est.method, " of 'shape'", sep = "")
    }, fourth.root = {
        p <- 1/4
        method <- paste("Normal Approximation of\n", space(33), 
            "Kulkarni & Powar (2010) using\n", space(33), "Hawkins & Wixley (1986) fourth-root\n", 
            space(33), "transformation to Normality", space(33), 
            "and ", shape.est.method, " of 'shape'", sep = "")
    })
    names(p) <- NULL
    y <- x^p
    n <- length(y)
    dum.list <- ci.norm(muhat = mean(y), sdhat = sd(y), n = length(y), 
        ci.type = ci.type, alpha = 1 - conf.level)
    LCL <- dum.list$limits["LCL"]
    UCL <- dum.list$limits["UCL"]
    if (ci.type %in% c("two-sided", "lower") & LCL < 0) {
        LCL <- 0
        warning("Normal approximation not accurate for this case")
    }
    switch(ci.type, `two-sided` = {
        LCL <- (LCL * shape^p * gamma(shape)/gamma(shape + p))^(1/p)
        UCL <- (UCL * shape^p * gamma(shape)/gamma(shape + p))^(1/p)
    }, lower = {
        LCL <- (LCL * shape^p * gamma(shape)/gamma(shape + p))^(1/p)
        UCL <- Inf
    }, upper = {
        LCL <- 0
        UCL <- (UCL * shape^p * gamma(shape)/gamma(shape + p))^(1/p)
    })
    ci.limits <- c(LCL, UCL)
    names(ci.limits) <- c("LCL", "UCL")
    interval <- list(name = "Confidence", parameter = "mean", 
        limits = ci.limits, type = ci.type, method = method, 
        conf.level = conf.level, sample.size = length(x), normal.transform.power = p)
    oldClass(interval) <- "intervalEstimate"
    interval
}
