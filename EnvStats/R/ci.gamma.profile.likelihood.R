ci.gamma.profile.likelihood <-
function (x, shape.mle, scale.mle, conf.level, LCL.start, UCL.start) 
{
    n <- length(x)
    mean.mle <- shape.mle * scale.mle
    cv.mle <- 1/sqrt(shape.mle)
    sd.mle <- sqrt(shape.mle)/scale.mle
    loglik.at.mle <- loglikComplete(theta = c(mean = mean.mle, 
        cv = cv.mle), x = x, distribution = "gammaAlt")
    fcn <- function(CL, loglik.at.mle, mean.mle, cv.mle, x) {
        cv.mle.at.CL <- egammaAlt.cv.mle.at.fixed.mean(fixed.mean = CL, 
            mean.mle = mean.mle, cv.mle = cv.mle, x = x)
        (2 * (loglik.at.mle - loglikComplete(theta = c(CL, cv.mle.at.CL), 
            x = x, distribution = "gammaAlt")) - qchisq(conf.level, 
            df = 1))^2
    }
    LCL <- nlminb(start = LCL.start, objective = fcn, lower = .Machine$double.eps, 
        upper = mean.mle, loglik.at.mle = loglik.at.mle, mean.mle = mean.mle, 
        cv.mle = cv.mle, x = x)$par
    UCL <- nlminb(start = UCL.start, objective = fcn, lower = mean.mle, 
        loglik.at.mle = loglik.at.mle, mean.mle = mean.mle, cv.mle = cv.mle, 
        x = x)$par
    ci.limits <- c(LCL, UCL)
    names(ci.limits) <- c("LCL", "UCL")
    interval <- list(name = "Confidence", parameter = "mean", 
        limits = ci.limits, type = "two-sided", method = "Profile Likelihood", 
        conf.level = conf.level)
    oldClass(interval) <- "intervalEstimate"
    interval
}
