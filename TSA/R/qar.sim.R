`qar.sim` <-
function (const = 0, phi0 = 0, phi1 = 0.5, sigma = 1, n = 20, 
    init = 0) 
{
    x = rep(0, n)
    x[1] = init
    e = rnorm(n) * sigma
    for (i in 2:n) x[i] = const + phi0 * x[i - 1] + phi1 * x[i - 
        1]^2 + e[i]
    x
}

