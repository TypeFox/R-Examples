chenTTest.ci <-
function (muhat, sdhat, skewhat, n, alternative, conf.level, 
    p.value.type = c("z", "t", "Avg. of z and t"), paired = FALSE) 
{
    p.value.type <- match.arg(p.value.type)
    alpha <- 1 - conf.level
    fcn.to.min <- function(mu.weird, muhat.weird, sdhat.weird, 
        skewhat.weird, n.weird, alternative.weird, alpha.weird, 
        p.value.type.weird) {
        p.value <- chenTTest.sub(mu = mu.weird, muhat = muhat.weird, 
            sdhat = sdhat.weird, skewhat = skewhat.weird, n = n.weird, 
            alternative = alternative.weird)$p.value[p.value.type.weird]
        (alpha.weird - p.value)^2
    }
    string <- ifelse(paired, "The sample skew for the paired differences", 
        "The sample skew")
    switch(alternative, greater = {
        ci.type <- "Lower"
        ucl <- Inf
        if (skewhat <= 0) {
            warning(paste(string, "is less than or equal to 0.\n ", 
                "Chen's test is not appropriate for a\n ", "\"lower\" confidence interval.\n"))
            lcl <- NA
        } else {
            lcl <- nlminb(start = muhat, objective = fcn.to.min, 
                upper = muhat, muhat.weird = muhat, sdhat.weird = sdhat, 
                skewhat.weird = skewhat, n.weird = n, alternative.weird = alternative, 
                alpha.weird = alpha, p.value.type.weird = p.value.type, 
                control = list(step.max = sdhat/sqrt(n)))$par
        }
    }, less = {
        ci.type <- "Upper"
        lcl <- -Inf
        if (skewhat >= 0) {
            warning(paste(string, "is greater than or equal to 0.\n ", 
                "Chen's test is not appropriate for an\n ", "\"upper\" confidence interval.\n"))
            ucl <- NA
        } else {
            ucl <- nlminb(start = muhat, objective = fcn.to.min, 
                lower = muhat, muhat.weird = muhat, sdhat.weird = sdhat, 
                skewhat.weird = skewhat, n.weird = n, alternative.weird = alternative, 
                alpha.weird = alpha, p.value.type.weird = p.value.type, 
                control = list(step.max = sdhat/sqrt(n)))$par
        }
    })
    ci.limits <- c(lcl, ucl)
    names(ci.limits) <- c("LCL", "UCL")
    ret.obj <- list(name = "Confidence", parameter = ifelse(paired, 
        "mean of differences", "mean"), limits = ci.limits, type = ci.type, 
        method = paste("Based on", p.value.type), conf.level = conf.level, 
        sample.size = n, dof = n - 1)
    oldClass(ret.obj) <- "intervalEstimate"
    ret.obj
}
