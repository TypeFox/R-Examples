egammaAlt.cv.mle.at.fixed.mean <-
function (fixed.mean, mean.mle, cv.mle, x) 
{
    fcn <- function(cv.mle.at.fixed.mean, fixed.mean, x) {
        -loglikComplete(theta = c(mean = fixed.mean, cv = cv.mle.at.fixed.mean), 
            x = x, distribution = "gammaAlt")
    }
    names(fixed.mean) <- NULL
    names(mean.mle) <- NULL
    names(cv.mle) <- NULL
    ret.val <- nlminb(start = cv.mle * mean.mle/fixed.mean, objective = fcn, 
        lower = 1e-07, fixed.mean = fixed.mean, x = x)$par
    names(ret.val) <- NULL
    ret.val
}
