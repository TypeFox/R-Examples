fix0106 <-
function (x, alpha = 4, beta = 4, mu = 0, sigma = 1, eta = 0, 
    kappa = 1) 
{
    eps <- (x - mu)/sigma
    ((eps/(2 * beta + eps^2)) * dnorm(x, mean = eta, sd = kappa))
}
