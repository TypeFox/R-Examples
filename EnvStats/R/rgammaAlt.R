rgammaAlt <-
function (n, mean, cv = 1) 
{
    shape <- cv^-2
    scale <- mean/shape
    stats::rgamma(n = n, shape = shape, scale = scale)
}
