## last modified June 2002

weibullparinv <- function(shape, scale, loc = 0) 
{
    nu <- 1/shape
    if (nu < 1e-6) {
        mu <- scale * (1 + nu * digamma(1) + nu^2 * (digamma(1)^2 + 
            trigamma(1))/2)
        sigma <- scale^2 * nu^2 * trigamma(1)
    }
    else {
        mu <- loc + gamma(1 + (nu)) * scale
        sigma <- sqrt(gamma(1 + 2 * nu) - (gamma(1 + nu))^2) * 
            scale
    }
    data.frame(mu, sigma, loc)
}
