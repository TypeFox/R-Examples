ci.normal.approx <-
function (theta.hat, sd.theta.hat, n, df = NULL, ci.type, alpha, 
    lb = -Inf, ub = Inf, test.statistic = c("t", "z")) 
{
    test.statistic <- match.arg(test.statistic)
    if (test.statistic == "t" && is.null(df)) 
        stop("When test.statistic='t' you must supply the 'df' argument")
    if (ci.type == "two.sided" || ci.type == "two-sided") {
        hw <- switch(test.statistic, t = qt(1 - (alpha/2), df) * 
            sd.theta.hat, z = qnorm(1 - (alpha/2)) * sd.theta.hat)
        lcl <- max(lb, theta.hat - hw)
        ucl <- min(ub, theta.hat + hw)
    }
    else {
        hw <- switch(test.statistic, t = qt(1 - alpha, df) * 
            sd.theta.hat, z = qnorm(1 - alpha) * sd.theta.hat)
        if (ci.type == "lower") {
            lcl <- theta.hat - hw
            ucl <- ub
        }
        else {
            lcl <- lb
            ucl <- theta.hat + hw
        }
    }
    ci.limits <- c(lcl, ucl)
    names(ci.limits) <- c("LCL", "UCL")
    method <- ifelse(test.statistic == "t", paste("Normal Approximation\n", 
        space(33), "(t Distribution)", sep = ""), "Normal Approximation")
    ret.obj <- list(name = "Confidence", parameter = "theta", 
        limits = ci.limits, type = ifelse(ci.type == "two.sided", 
            "two-sided", ci.type), method = method, conf.level = 1 - 
            alpha)
    if (test.statistic == "t") 
        ret.obj <- c(ret.obj, list(sample.size = n, dof = df))
    oldClass(ret.obj) <- "intervalEstimate"
    ret.obj
}
