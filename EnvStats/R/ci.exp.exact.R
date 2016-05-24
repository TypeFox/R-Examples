ci.exp.exact <-
function (rate, n, ci.type, alpha) 
{
    df <- 2 * n
    denom <- df/rate
    switch(ci.type, `two-sided` = {
        lcl <- qchisq(alpha/2, df)/denom
        ucl <- qchisq(1 - alpha/2, df)/denom
    }, lower = {
        lcl <- qchisq(alpha, df)/denom
        ucl <- Inf
    }, upper = {
        lcl <- 0
        ucl <- qchisq(1 - alpha, df)/denom
    })
    ci.limits <- c(lcl, ucl)
    names(ci.limits) <- c("LCL", "UCL")
    ret.obj <- list(name = "Confidence", parameter = "rate", 
        limits = ci.limits, type = ci.type, method = "Exact", 
        conf.level = 1 - alpha, sample.size = n)
    oldClass(ret.obj) <- "intervalEstimate"
    ret.obj
}
