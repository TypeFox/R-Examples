dgammaAlt <-
function (x, mean, cv = 1, log = FALSE) 
{
    shape <- cv^-2
    scale <- mean/shape
    stats::dgamma(x = x, shape = shape, scale = scale, log = log)
}
