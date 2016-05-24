fix0108 <-
function (func = 1, a = 0, b = 1, n = 1000, alpha, beta, mu, 
    sigma, eta, kappa) 
{
    (sum(func(a + ((1:n - 0.5) * (b - a))/n, alpha, beta, mu, 
        sigma, eta, kappa)) * (b - a))/n
}
