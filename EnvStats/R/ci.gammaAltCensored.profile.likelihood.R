ci.gammaAltCensored.profile.likelihood <-
function (x, censored, censoring.side, mean.mle, cv.mle, ci.type, 
    conf.level, LCL.normal.approx, UCL.normal.approx) 
{
    loglik.at.mle <- loglikCensored(theta = c(mean.mle, cv.mle), 
        x = x, censored = censored, censoring.side = censoring.side, 
        distribution = "gammaAlt")
    fcn <- function(CL, loglik.at.mle, cv.mle, x, censored, censoring.side, 
        conf.level) {
        cv.mle.at.CL <- egammaAltCensored.cv.mle.at.fixed.mean(fixed.mean = CL, 
            mean.mle = mean.mle, cv.mle = cv.mle, x = x, censored = censored, 
            censoring.side = censoring.side)
        (2 * (loglik.at.mle - loglikCensored(theta = c(CL, cv.mle.at.CL), 
            x = x, censored = censored, censoring.side = censoring.side, 
            distribution = "gammaAlt")) - qchisq(conf.level, 
            df = 1))^2
    }
    switch(ci.type, `two-sided` = {
        LCL <- nlminb(start = LCL.normal.approx, objective = fcn, 
            lower = .Machine$double.eps, upper = mean.mle, loglik.at.mle = loglik.at.mle, 
            cv.mle = cv.mle, x = x, censored = censored, censoring.side = censoring.side, 
            conf.level = conf.level)$par
        UCL <- nlminb(start = UCL.normal.approx, objective = fcn, 
            lower = mean.mle, loglik.at.mle = loglik.at.mle, 
            cv.mle = cv.mle, x = x, censored = censored, censoring.side = censoring.side, 
            conf.level = conf.level)$par
    }, lower = {
        LCL <- nlminb(start = LCL.normal.approx, objective = fcn, 
            lower = .Machine$double.eps, upper = mean.mle, loglik.at.mle = loglik.at.mle, 
            cv.mle = cv.mle, x = x, censored = censored, censoring.side = censoring.side, 
            conf.level = 1 - 2 * (1 - conf.level))$par
        UCL <- Inf
    }, upper = {
        LCL <- 0
        UCL <- nlminb(start = UCL.normal.approx, objective = fcn, 
            lower = mean.mle, loglik.at.mle = loglik.at.mle, 
            cv.mle = cv.mle, x = x, censored = censored, censoring.side = censoring.side, 
            conf.level = 1 - 2 * (1 - conf.level))$par
    })
    ci.limits <- c(LCL, UCL)
    names(ci.limits) <- c("LCL", "UCL")
    interval <- list(name = "Confidence", parameter = "mean", 
        limits = ci.limits, type = ci.type, method = "Profile Likelihood", 
        conf.level = conf.level)
    oldClass(interval) <- "intervalEstimate"
    interval
}
