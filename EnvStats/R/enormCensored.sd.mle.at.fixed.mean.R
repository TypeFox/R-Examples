enormCensored.sd.mle.at.fixed.mean <-
function (fixed.mean, sd.mle, x, censored, censoring.side) 
{
    fcn <- function(sd.mle.at.fixed.mean, fixed.mean, x, censored, 
        censoring.side) {
        -loglikCensored(theta = c(mean = fixed.mean, sd = sd.mle.at.fixed.mean), 
            x = x, censored = censored, censoring.side = censoring.side, 
            distribution = "norm")
    }
    names(fixed.mean) <- NULL
    names(sd.mle) <- NULL
    ret.val <- nlminb(start = sd.mle, objective = fcn, lower = .Machine$double.eps, 
        fixed.mean = fixed.mean, x = x, censored = censored, 
        censoring.side = censoring.side)$par
    names(ret.val) <- NULL
    ret.val
}
