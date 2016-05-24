ci.lnorm3.zero.skew <-
function (threshold, meanlog, sdlog, x, x1, n, threshold.lb, 
    threshold.ub, ci.type, alpha, ci.parameter) 
{
    conf.level <- 1 - alpha
    df <- n - 2
    fcn.to.optimize <- function(threshold, x.weird, alternative, 
        conf.level) {
        (skewGofTest(log(x.weird - threshold), alternative = alternative)$p.value - 
            (1 - conf.level))^2
    }
    switch(ci.type, `two-sided` = {
        nlminb.list.lcl <- nlminb(start = max(threshold - 1, 
            threshold.lb), objective = fcn.to.optimize, lower = threshold.lb, 
            upper = threshold - .Machine$double.eps, x.weird = x, 
            alternative = "greater", conf.level = 1 - alpha/2)
        lcl <- nlminb.list.lcl$par
        nlminb.list.ucl <- nlminb(start = mean(threshold, x1), 
            objective = fcn.to.optimize, lower = threshold + 
                .Machine$double.eps, upper = threshold.ub, x.weird = x, 
            alternative = "less", conf.level = 1 - alpha/2)
        ucl <- nlminb.list.ucl$par
        if (nlminb.list.lcl$objective > .Machine$double.eps || 
            nlminb.list.ucl$objective > .Machine$double.eps || 
            lcl == threshold.lb || ucl == threshold.ub) {
            warning(paste("Unable to find confidence limits for 'threshold'", 
                "in the interval", "[ mean(x) - threshold.lb.sd * sd(x), min(x) ). ", 
                "Try changing the value of 'threshold.lb.sd'"))
            lcl <- ucl <- NA
        }
    }, lower = {
        nlminb.list.lcl <- nlminb(start = max(threshold - 1, 
            threshold.lb), objective = fcn.to.optimize, lower = threshold.lb, 
            upper = threshold - .Machine$double.eps, x.weird = x, 
            alternative = "greater", conf.level = conf.level)
        lcl <- nlminb.list.lcl$par
        ucl <- x1
        if (nlminb.list.lcl$objective > .Machine$double.eps || 
            lcl == threshold.lb) {
            warning(paste("Unable to find lower confidence limit", 
                "for 'threshold' in the interval", "[ mean(x) - 'threshold.lb.sd' * sd(x), 'threshold.hat' ).", 
                "Try changing the value of 'threshold.lb.sd'"))
            lcl <- NA
        }
    }, upper = {
        lcl <- -Inf
        nlminb.list.ucl <- nlminb(start = mean(threshold, x1), 
            objective = fcn.to.optimize, lower = threshold + 
                .Machine$double.eps, upper = threshold.ub, x.weird = x, 
            alternative = "less", conf.level = conf.level)
        ucl <- nlminb.list.ucl$par
        if (nlminb.list.ucl$objective > .Machine$double.eps || 
            ucl == threshold.ub) {
            warning(paste("Unable to find upper confidence limit", 
                "for 'threshold' in the interval", "('threshold.hat' , min(x) )"))
            ucl <- NA
        }
    })
    ci.limits <- c(lcl, ucl)
    names(ci.limits) <- c("LCL", "UCL")
    if (ci.parameter == "median") {
        ci.obj.meanlog <- ci.normal.approx(meanlog, sdlog/sqrt(n), 
            n = n, df = n - 2, ci.type = ci.type, alpha = alpha)
        ci.limits <- ci.limits + exp(ci.obj.meanlog$limits)
    }
    ret.obj <- list(name = "Confidence", parameter = ci.parameter, 
        limits = ci.limits, type = ci.type, method = "Skewness", 
        conf.level = conf.level, sample.size = n, dof = df)
    oldClass(ret.obj) <- "intervalEstimate"
    ret.obj
}
