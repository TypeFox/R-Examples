ci.qpois <-
function (p, lambda.hat, conf.limits, n, ci.type, alpha, digits = 0) 
{
    conf.level <- 1 - alpha
    switch(ci.type, `two-sided` = {
        lcl <- qpois(p, conf.limits["LCL"])
        ucl <- qpois(p, conf.limits["UCL"])
    }, lower = {
        lcl <- qpois(p, conf.limits["LCL"])
        ucl <- Inf
    }, upper = {
        lcl <- 0
        ucl <- qpois(p, conf.limits["UCL"])
    })
    ci.limits <- c(lcl, ucl)
    names(ci.limits) <- c("LCL", "UCL")
    ret.obj <- list(name = "Confidence", parameter = paste(round(100 * 
        p, digits), "'th %ile", sep = ""), limits = ci.limits, 
        type = ci.type, method = "Exact", conf.level = conf.level, 
        sample.size = n)
    oldClass(ret.obj) <- "intervalEstimate"
    ret.obj
}
