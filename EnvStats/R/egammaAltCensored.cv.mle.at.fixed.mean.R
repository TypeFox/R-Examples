egammaAltCensored.cv.mle.at.fixed.mean <-
function (fixed.mean, mean.mle, cv.mle, x, censored, censoring.side) 
{
    fcn <- function(cv.mle.at.fixed.mean, fixed.mean, x, censored, 
        censoring.side) {
        -loglikCensored(theta = c(mean = fixed.mean, cv = cv.mle.at.fixed.mean), 
            x = x, censored = censored, censoring.side = censoring.side, 
            distribution = "gammaAlt")
    }
    names(fixed.mean) <- NULL
    names(mean.mle) <- NULL
    names(cv.mle) <- NULL
    ret.val <- nlminb(start = cv.mle * mean.mle/fixed.mean, objective = fcn, 
        lower = 1e-07, fixed.mean = fixed.mean, x = x, censored = censored, 
        censoring.side = censoring.side)$par
    names(ret.val) <- NULL
    ret.val
}
