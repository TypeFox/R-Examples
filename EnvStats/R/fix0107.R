fix0107 <-
function (func = 1, a = 0, b = 1, n = 1000, alpha, beta, mu, 
    sigma, eta, kappa) 
{
    i0 <- fix0108(func, a, b, n, alpha, beta, mu, sigma, eta, 
        kappa)
    for (i in (2:100) * 1000) {
        i1 <- fix0108(func, a, b, i, alpha, beta, mu, sigma, 
            eta, kappa)
        if (abs(i1 - i0)/i1 < 1e-05) 
            break
        i0 <- i1
    }
    i1
}
