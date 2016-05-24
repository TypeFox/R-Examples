ci.norm.var <-
function (sdhat, n, ci.type, alpha) 
{
    varhat <- sdhat^2
    df <- n - 1
    if (ci.type == "two.sided" || ci.type == "two-sided") {
        lcl <- (df * varhat)/qchisq(1 - alpha/2, df)
        ucl <- (df * varhat)/qchisq(alpha/2, df)
    }
    else if (ci.type == "lower") {
        lcl <- (df * varhat)/qchisq(1 - alpha, df)
        ucl <- Inf
    }
    else {
        lcl <- 0
        ucl <- (df * varhat)/qchisq(alpha, df)
    }
    ci.limits <- c(lcl, ucl)
    names(ci.limits) <- c("LCL", "UCL")
    method <- "Exact"
    ret.obj <- list(name = "Confidence", parameter = "variance", 
        limits = ci.limits, type = ifelse(ci.type == "two.sided", 
            "two-sided", ci.type), method = method, conf.level = 1 - 
            alpha, sample.size = n, dof = df)
    oldClass(ret.obj) <- "intervalEstimate"
    ret.obj
}
