.rinvgamma <- function(n, shape = 1, scale = 1)
{
    return(1 / stats::rgamma(n, shape, rate = scale))
}
